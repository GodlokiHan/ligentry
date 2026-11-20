"""Dual-tier PyTorch pose-rescoring architectures (CPU-light & GPU-deep).

Both models accept identical fixed-shape padded bipartite graph tensors so a single
ONNX Runtime plugin in Ligentry can load either `cpu_v1.onnx` or `gpu_v1.onnx`.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F

from ligentry.plugins._rescore_features import EDGE_FEAT_DIM, NODE_FEAT_DIM


class EdgeMPNNBlock(nn.Module):
    """Edge-conditioned message passing block with residual connection & LayerNorm."""

    def __init__(self, hidden_dim: int, edge_dim: int, dropout: float = 0.1) -> None:
        super().__init__()
        self.msg_mlp = nn.Sequential(
            nn.Linear(hidden_dim * 2 + edge_dim, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.gate = nn.Sequential(
            nn.Linear(edge_dim, hidden_dim),
            nn.Sigmoid(),
        )
        self.update_mlp = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.norm = nn.LayerNorm(hidden_dim)

    def forward(
        self,
        h: torch.Tensor,          # [B, N, D]
        edge_index: torch.Tensor, # [B, 2, E]
        edge_attr: torch.Tensor,  # [B, E, De]
        edge_mask: torch.Tensor,  # [B, E]
    ) -> torch.Tensor:
        bsz, n_nodes, dim = h.shape
        src = edge_index[:, 0, :]  # [B, E]
        dst = edge_index[:, 1, :]  # [B, E]

        batch_idx = torch.arange(bsz, device=h.device)[:, None]
        h_src = h[batch_idx, src]  # [B, E, D]
        h_dst = h[batch_idx, dst]  # [B, E, D]

        msg_in = torch.cat([h_src, h_dst, edge_attr], dim=-1)
        msg = self.msg_mlp(msg_in) * self.gate(edge_attr)
        msg = msg * edge_mask.unsqueeze(-1)

        agg = torch.zeros_like(h)
        dst_exp = dst.unsqueeze(-1).expand(-1, -1, dim)
        agg.scatter_add_(1, dst_exp, msg)

        out = self.update_mlp(torch.cat([h, agg], dim=-1))
        return self.norm(h + out)


class LigandAttentionPool(nn.Module):
    """Attention-weighted pooling restricted to ligand heavy atoms."""

    def __init__(self, hidden_dim: int, heads: int = 4) -> None:
        super().__init__()
        self.attn = nn.Linear(hidden_dim, heads)
        self.proj = nn.Linear(hidden_dim * heads, hidden_dim)

    def forward(
        self,
        h: torch.Tensor,         # [B, N, D]
        node_mask: torch.Tensor, # [B, N]
        is_ligand: torch.Tensor, # [B, N]
    ) -> torch.Tensor:
        lig_mask = (node_mask > 0.5) & (is_ligand > 0.5)  # [B, N]
        logits = self.attn(h)  # [B, N, H]
        logits = logits.masked_fill(~lig_mask.unsqueeze(-1), torch.finfo(logits.dtype).min)
        weights = torch.softmax(logits, dim=1)  # [B, N, H]
        pooled = torch.bmm(weights.transpose(1, 2), h).reshape(h.size(0), -1)
        return self.proj(pooled)


class LigentryRescorerNet(nn.Module):
    """Unified pose rescoring network parameterized for CPU-default or GPU-heavy tiers."""

    def __init__(
        self,
        hidden_dim: int = 96,
        num_layers: int = 3,
        pool_heads: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.node_in = nn.Sequential(
            nn.Linear(NODE_FEAT_DIM, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
        )
        self.edge_in = nn.Sequential(
            nn.Linear(EDGE_FEAT_DIM, hidden_dim // 2),
            nn.SiLU(),
            nn.Linear(hidden_dim // 2, hidden_dim // 2),
        )
        self.blocks = nn.ModuleList(
            [EdgeMPNNBlock(hidden_dim, hidden_dim // 2, dropout=dropout) for _ in range(num_layers)]
        )
        self.lig_pool = LigandAttentionPool(hidden_dim, heads=pool_heads)

        self.head = nn.Sequential(
            nn.Linear(hidden_dim + 1, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 2),
        )
        self.vina_skip = nn.Parameter(torch.tensor([-0.25]))

    def forward(
        self,
        node_x: torch.Tensor,     # [P, N, 16]
        node_mask: torch.Tensor,  # [P, N]
        edge_index: torch.Tensor, # [P, 2, E]
        edge_attr: torch.Tensor,  # [P, E, 20]
        edge_mask: torch.Tensor,  # [P, E]
        vina_score: torch.Tensor, # [P, 1]
    ) -> tuple[torch.Tensor, torch.Tensor]:
        is_ligand = node_x[:, :, 15]
        h = self.node_in(node_x) * node_mask.unsqueeze(-1)
        e = self.edge_in(edge_attr) * edge_mask.unsqueeze(-1)

        for block in self.blocks:
            h = block(h, edge_index, e, edge_mask)
            h = h * node_mask.unsqueeze(-1)

        g = self.lig_pool(h, node_mask, is_ligand)
        v_norm = torch.clamp(vina_score / 10.0, -2.5, 2.5)
        out = self.head(torch.cat([g, v_norm], dim=-1))

        pose_logit = out[:, 0:1] + self.vina_skip * vina_score
        pkd_pred = F.softplus(out[:, 1:2])
        return pose_logit, pkd_pred


def build_model(tier: str = "cpu") -> LigentryRescorerNet:
    if tier == "cpu":
        return LigentryRescorerNet(hidden_dim=96, num_layers=3, pool_heads=4, dropout=0.10)
    if tier == "gpu":
        return LigentryRescorerNet(hidden_dim=192, num_layers=5, pool_heads=8, dropout=0.15)
    raise ValueError(f"Unknown tier '{tier}', expected 'cpu' or 'gpu'")

