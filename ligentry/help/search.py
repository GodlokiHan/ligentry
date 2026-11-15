"""Search the manual (no Qt).

Every word you type must appear somewhere in a topic (title, keywords, summary or
body, prefix matches allowed) for the topic to be listed; topics are ranked by where
the words were found, with the title counting most.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from ligentry.help.topics import TOPICS, Topic

_TAG = re.compile(r"<[^>]+>")
_WORD = re.compile(r"[a-z0-9][a-z0-9+\-./]*", re.IGNORECASE)
_ENTITIES = {"&amp;": "&", "&lt;": "<", "&gt;": ">", "&quot;": '"', "&Aring;": "A", "&times;": "x",
             "&middot;": " ", "&plusmn;": "+-", "&minus;": "-", "&nbsp;": " "}


def plain_text(html: str) -> str:
    """The text of an HTML snippet, with tags and common entities removed."""
    text = _TAG.sub(" ", html)
    for ent, rep in _ENTITIES.items():
        text = text.replace(ent, rep)
    return " ".join(text.split())


def tokens(text: str) -> list[str]:
    return [w.strip(".-/").lower() for w in _WORD.findall(text) if w.strip(".-/")]


@dataclass(frozen=True)
class Hit:
    topic: Topic
    score: int
    snippet: str


def _fields(topic: Topic) -> dict[str, str]:
    return {
        "title": topic.title.lower(),
        "keywords": " ".join(topic.keywords).lower(),
        "summary": topic.summary.lower(),
        "body": plain_text(topic.body).lower(),
    }


_WEIGHT = {"title": 10, "keywords": 6, "summary": 4, "body": 1}


def _word_hits(field_text: str, word: str) -> int:
    """How many words in the text start with `word`."""
    return sum(1 for w in tokens(field_text) if w.startswith(word))


def _snippet(topic: Topic, words: list[str], width: int = 140) -> str:
    text = plain_text(topic.body)
    low = text.lower()
    for w in words:
        i = low.find(w)
        if i >= 0:
            start = max(0, i - width // 3)
            piece = text[start:start + width].strip()
            return ("..." if start > 0 else "") + piece + ("..." if start + width < len(text) else "")
    return topic.summary


def search(query: str, limit: int = 20) -> list[Hit]:
    """Topics matching every word of `query`, best first. An empty query gives no hits."""
    words = tokens(query)
    if not words:
        return []
    hits: list[Hit] = []
    for topic in TOPICS:
        fields = _fields(topic)
        score = 0
        ok = True
        for w in words:
            word_score = 0
            for name, text in fields.items():
                n = _word_hits(text, w)
                if n:
                    word_score += _WEIGHT[name] * (1 + min(n, 5) * 0.1) if name == "body" else _WEIGHT[name]
            if word_score == 0:
                ok = False
                break
            score += word_score
        if ok:
            hits.append(Hit(topic, int(round(score)), _snippet(topic, words)))
    hits.sort(key=lambda h: (-h.score, h.topic.title))
    return hits[:limit]
