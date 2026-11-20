"""Tests for the manual, search, web helpers and welcome guide (no Qt needed).

Run from the project folder:  python scripts/test_help.py
"""
from __future__ import annotations

import re
import sys
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ligentry.help import guide, pages, search, topics, web  # noqa: E402

passed = failed = 0


def check(name: str, cond: bool, detail: str = "") -> None:
    global passed, failed
    if cond:
        passed += 1
    else:
        failed += 1
        print(f"FAIL: {name} {detail}")


def raises(fn, exc=ValueError) -> bool:
    try:
        fn()
    except exc:
        return True
    except Exception:
        return False
    return False


# ---------------------------------------------------------------- topics: structure
ids = [t.id for t in topics.TOPICS]
check("topic ids unique", len(ids) == len(set(ids)))
check("home id not used by a topic", pages.HOME_ID not in ids)
check("at least 30 topics", len(ids) >= 30, str(len(ids)))
for t in topics.TOPICS:
    check(f"{t.id}: section known", t.section in topics.SECTIONS, t.section)
    check(f"{t.id}: title/summary/body present", bool(t.title and t.summary and t.body.strip()))
    check(f"{t.id}: summary is one line", "\n" not in t.summary and "<" not in t.summary)
    for r in t.related:
        check(f"{t.id}: related '{r}' exists", topics.has_topic(r))
        check(f"{t.id}: does not list itself", r != t.id)
    check(f"{t.id}: web queries are searchable", all(web.build_search_url("web", q) for q in t.web))
for s in topics.SECTIONS:
    check(f"section '{s}' not empty", len(topics.topics_in_section(s)) > 0)
check("get_topic unknown raises", raises(lambda: topics.get_topic("nope"), KeyError))

# ---------------------------------------------------------------- topics: HTML safety
ALLOWED = {"h2", "p", "ul", "ol", "li", "b", "i", "code", "kbd", "pre", "h3", "h4", "table", "tr", "th", "td", "br", "a", "sup"}
VOID = {"br"}


class Checker(HTMLParser):
    def __init__(self):
        super().__init__()
        self.stack: list[str] = []
        self.errors: list[str] = []
        self.links: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag not in ALLOWED:
            self.errors.append(f"tag <{tag}> not allowed")
        for k, v in attrs:
            if k == "style" or k == "color":
                self.errors.append("style/color attribute")
            if tag == "a" and k == "href":
                self.links.append(v or "")
        if tag not in VOID:
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if not self.stack or self.stack[-1] != tag:
            self.errors.append(f"unbalanced </{tag}> (open: {self.stack[-3:]})")
        else:
            self.stack.pop()

    def close(self):
        super().close()
        if self.stack:
            self.errors.append(f"unclosed {self.stack}")


all_pages = {t.id: pages.render_page(t.id) for t in topics.TOPICS}
all_pages[pages.HOME_ID] = pages.render_home()
for pid, body in all_pages.items():
    c = Checker()
    c.feed(body)
    c.close()
    check(f"{pid}: html clean", not c.errors, "; ".join(c.errors[:3]))
    for href in c.links:
        link = pages.classify_link(href)
        check(f"{pid}: link {href[:50]} resolves", link.kind in ("topic", "web"), f"{link.kind}")
        if link.kind == "web":
            check(f"{pid}: web link safe {href[:40]}", web.is_safe_url(link.value))
for g in guide.GUIDE_PAGES:
    c = Checker()
    c.feed(g.body)
    c.close()
    check(f"guide '{g.title}': html clean", not c.errors, "; ".join(c.errors[:3]))
    check(f"guide '{g.title}': topic exists", topics.has_topic(g.topic))

# ---------------------------------------------------------------- every topic reachable
linked = set()
for t in topics.TOPICS:
    for m in re.finditer(r'href="topic:([a-z0-9-]+)"', t.body):
        linked.add(m.group(1))
    linked.update(t.related)
home = pages.render_home()
check("every topic is on the contents page", all(f'topic:{t.id}"' in home for t in topics.TOPICS))

# ---------------------------------------------------------------- web helpers
check("source keys unique", len({s.key for s in web.SOURCES}) == len(web.SOURCES))
for s in web.SOURCES:
    check(f"{s.key}: template has {{q}} once", s.template.count("{q}") == 1)
    url = web.build_search_url(s.key, "ke0 effect compartment")
    check(f"{s.key}: url safe", web.is_safe_url(url), url)
    check(f"{s.key}: https", url.startswith("https://"))
    check(f"{s.key}: query encoded", " " not in url)
check("default source exists", web.get_source(web.DEFAULT_SOURCE) is not None)
check("quick sources exist", all(web.get_source(k) for k in web.QUICK_SOURCES))
check("unknown source", raises(lambda: web.get_source("zzz")))
check("empty query", raises(lambda: web.build_search_url("pubmed", "   ")))
check("long query", raises(lambda: web.build_search_url("pubmed", "x" * 501)))
check("special chars encoded", "%26" in web.build_search_url("pubmed", "a&b") and "&b" not in web.build_search_url("pubmed", "a&b"))
check("hash chars encoded", "%23" in web.build_search_url("pubchem", "a#b"))
check("whitespace collapsed", web.build_search_url("pubmed", "a   b\n c").endswith("a+b+c"))

for url, ok in [
    ("https://pubmed.ncbi.nlm.nih.gov/123/", True),
    ("http://example.org", True),
    ("https://8.8.8.8/x", True),
    ("file:///etc/passwd", False),
    ("javascript:alert(1)", False),
    ("ftp://example.org", False),
    ("data:text/html,hi", False),
    ("https://localhost/x", False),
    ("http://127.0.0.1:8000", False),
    ("http://192.168.1.1", False),
    ("http://10.0.0.5", False),
    ("http://169.254.169.254/latest", False),
    ("http://[::1]/", False),
    ("http://printer.local/", False),
    ("http://intranet/", False),
    ("https://user:pw@example.org/", False),
    ("https://", False),
    ("", False),
]:
    check(f"is_safe_url {url!r}", web.is_safe_url(url) == ok)

R = web.resolve_input
check("doi bare", R("10.1000/xyz123").url == "https://doi.org/10.1000/xyz123")
check("doi prefixed", R("doi: 10.1000/xyz123").url == "https://doi.org/10.1000/xyz123")
check("doi url", R("https://doi.org/10.1000/xyz123").url == "https://doi.org/10.1000/xyz123")
check("doi trailing period", R("10.1000/xyz123.").url == "https://doi.org/10.1000/xyz123")
check("pmid", R("PMID 12345678").url == "https://pubmed.ncbi.nlm.nih.gov/12345678/")
check("pmid colon", R("pmid:42").url == "https://pubmed.ncbi.nlm.nih.gov/42/")
check("pmc", R("pmc1234567").url == "https://pmc.ncbi.nlm.nih.gov/articles/PMC1234567/")
check("nct", R("nct01234567").url == "https://clinicaltrials.gov/study/NCT01234567")
check("pdb", R("PDB 1hsg").url == "https://www.rcsb.org/structure/1HSG")
check("domain gets https", R("who.int").url == "https://who.int")
check("domain with path", R("pubmed.ncbi.nlm.nih.gov/123").url == "https://pubmed.ncbi.nlm.nih.gov/123")
check("full url kept", R("https://www.fda.gov/drugs").url == "https://www.fda.gov/drugs")
check("localhost refused", raises(lambda: R("http://localhost:8080")))
check("private ip refused", raises(lambda: R("http://192.168.0.1/admin")))
check("empty refused", raises(lambda: R("  ")))
for text in ["1.5", "e.g.", "ke0", "half life of warfarin", "file:///etc/passwd", "javascript:alert(1)"]:
    r = R(text, "pubmed")
    check(f"plain text is a search: {text!r}", r.is_search and r.url.startswith("https://pubmed.ncbi.nlm.nih.gov/?term="), r.url)
check("search uses chosen source", R("aspirin", "youtube").url.startswith("https://www.youtube.com/results"))
check("search default source", R("aspirin").url.startswith("https://duckduckgo.com/"))
check("search label", "PubMed" in R("aspirin", "pubmed").label)

# ---------------------------------------------------------------- manual search
check("empty query no hits", search.search("") == [] and search.search("   ") == [])
check("punctuation only no hits", search.search("?!") == [])
check("no match", search.search("zzzqqq") == [])
first = lambda q: search.search(q)[0].topic.id if search.search(q) else None  # noqa: E731
check("search ke0", first("ke0") in ("pkpd-delayed", "pkpd-effect-fit"), str(first("ke0")))
check("search 'delayed effect'", first("delayed effect") == "pkpd-delayed", str(first("delayed effect")))
check("search stereoisomers", first("stereoisomers") == "docking-stereo", str(first("stereoisomers")))
check("search 'half-life'", search.search("half-life") != [])
check("search shortcuts", first("keyboard shortcuts") == "shortcuts", str(first("keyboard shortcuts")))
check("search benchmark", first("benchmark") == "benchmark", str(first("benchmark")))
check("search smiles", "docking-inputs" in [h.topic.id for h in search.search("smiles")])
check("search flip-flop", "pkpd-fit" in [h.topic.id for h in search.search("flip-flop")])
check("prefix match", first("stereo") in ("docking-stereo",), str(first("stereo")))
check("case insensitive", [h.topic.id for h in search.search("KE0")] == [h.topic.id for h in search.search("ke0")])
check("all words must match", search.search("ke0 zzzqqq") == [])
check("limit", len(search.search("the", limit=3)) <= 3)
check("hits sorted", all(a.score >= b.score for a, b in zip(search.search("model"), search.search("model")[1:])))
check("snippet not empty", all(h.snippet for h in search.search("effect")))
check("snippet is plain text", all("<" not in h.snippet for h in search.search("effect")))
check("plain_text entities", search.plain_text("<p>a &lt; b &amp; c</p>") == "a < b & c")
# every topic can be found by its own title
for t in topics.TOPICS:
    found = [h.topic.id for h in search.search(t.title)]
    check(f"{t.id}: found by title", t.id in found[:5], str(found[:5]))

# ---------------------------------------------------------------- pages and history
check("search page for no hits", "Nothing in the manual" in pages.render_search_results("zzz", []))
check("search page escapes input", "<script>" not in pages.render_search_results("<script>", []))
check("search page lists hits", "docking-stereo" in pages.render_search_results("x", search.search("stereoisomers")))
check("render home", pages.render("home").startswith("<h2>ligentry manual"))
check("render topic", "<h2>Stereoisomers</h2>" in pages.render("docking-stereo"))
check("render search", "Search results" in pages.render("search:ke0", search.search("ke0")))
check("topic page has web section", "Look it up on the web" in pages.render_page("pkpd-delayed"))
check("topic page has see also", "See also" in pages.render_page("pkpd-delayed"))
check("web link classified", pages.classify_link(f"websearch:pubmed:{'a%20b'}").value == "https://pubmed.ncbi.nlm.nih.gov/?term=a+b")
check("topic link classified", pages.classify_link("topic:glossary").kind == "topic")
check("home link classified", pages.classify_link("topic:home").kind == "topic")
check("bad topic link", pages.classify_link("topic:nope").kind == "unknown")
check("bad websearch", pages.classify_link("websearch:zzz:abc").kind == "unknown")
check("file link blocked/unknown", pages.classify_link("file:///etc/passwd").kind == "unknown")
check("localhost link blocked", pages.classify_link("http://localhost/x").kind == "blocked")
check("https link", pages.classify_link("https://www.rcsb.org/").kind == "web")

h = pages.History()
check("history start", h.current == "home" and not h.can_back and not h.can_forward)
h.visit("glossary")
h.visit("units")
check("history visit", h.current == "units" and h.can_back and not h.can_forward)
check("history back", h.back() == "glossary" and h.can_forward)
check("history forward", h.forward() == "units")
h.back()
h.visit("cli")
check("visit drops forward", h.current == "cli" and not h.can_forward and h.back() == "glossary")
before = h.current
h.visit(before)
check("revisit same page no-op", h.back() == "home")
check("back at start stays", h.back() == "home")

# ---------------------------------------------------------------- guide logic


class FakeSettings(dict):
    def value(self, key, default=None):
        return self.get(key, default)

    def setValue(self, key, value):
        self[key] = value


s = FakeSettings()
check("guide shown on first launch", guide.should_show_guide(s, {}))
guide.mark_guide_seen(s)
check("guide not shown again", not guide.should_show_guide(s, {}))
guide.mark_guide_seen(s, show_again=True)
check("show again honoured", guide.should_show_guide(s, {}))
guide.mark_guide_seen(s)
s[guide.SETTINGS_KEY] = "garbage"
check("corrupt setting shows guide", guide.should_show_guide(s, {}))
s2 = FakeSettings()
check("env suppresses", not guide.should_show_guide(s2, {guide.ENV_SUPPRESS: "1"}))
s3 = FakeSettings({guide.SETTINGS_KEY: str(guide.GUIDE_VERSION)})
check("string setting value ok (QSettings returns str)", not guide.should_show_guide(s3, {}))
check("newer guide version shows again", guide.should_show_guide(FakeSettings({guide.SETTINGS_KEY: guide.GUIDE_VERSION - 1}), {}))
check("guide has pages", len(guide.GUIDE_PAGES) >= 3)

# ---------------------------------------------------------------- drift: manual vs program
APP = ROOT / "ligentry" / "gui" / "app.py"
PK = ROOT / "ligentry" / "gui" / "pkpd_window.py"
HM = ROOT / "ligentry" / "gui" / "help_menu.py"
if APP.exists():
    app_src = APP.read_text(encoding="utf-8")
    shortcuts_body = topics.get_topic("shortcuts").body
    for key in ["Ctrl+N", "Ctrl+O", "Ctrl+S", "Ctrl+Shift+S", "Ctrl+E", "Ctrl+Q", "F2", "Ctrl+Shift+C",
                "Ctrl+Backspace", "Ctrl+R", "Ctrl+.", "Ctrl+Shift+P"]:
        check(f"app.py defines {key}", f'"{key}"' in app_src)
        check(f"manual lists {key}", f"<td>{key}</td>" in shortcuts_body)
    for pat, label in [(r"n_spin.*?setValue\(5\)", "max ligands 5"), (r"exh_spin.*?setValue\(8\)", "exhaustiveness 8"),
                       (r"iso_spin.*?setValue\(8\)", "max isomers 8")]:
        check(f"app.py default {label}", re.search(pat, app_src, re.S) is not None)
    check("app.py default PDB 1HSG", "1HSG" in app_src and "CHEMBL243" in app_src)
    check("app.py spin ranges match manual", "setRange(1, 100_000)" in app_src and "setRange(1, 64)" in app_src)
    check("manual says runs/gui", "runs/gui" in topics.get_topic("docking-run").body and "runs/gui" in app_src)
else:
    print("note: ligentry/gui/app.py not found, skipping drift checks against it")
if HM.exists():
    hm_src = HM.read_text(encoding="utf-8")
    for key in ["F1", "Ctrl+Shift+F", "Ctrl+L", "Alt+Left", "Alt+Right", "Alt+Home", "Ctrl+W"]:
        check(f"help code defines {key}", key in hm_src or key in (ROOT / "ligentry" / "gui" / "help_window.py").read_text(encoding="utf-8"))
if PK.exists():
    pk_src = PK.read_text(encoding="utf-8")
    check("pkpd window has export CSV text", "Export simulation CSV" in pk_src)
    check("pkpd window has Delayed effect text", "Delayed effect" in pk_src)
    check("pkpd window has Fit effect", "Fit effect" in pk_src and "Use fitted effect values" in pk_src)

# CSV columns quoted in the manual match the code
import importlib  # noqa: E402

try:
    runio_cols = None
    wf = importlib.import_module("ligentry.pkpd.workflow")
    src = Path(wf.__file__).read_text(encoding="utf-8")
    export_body = topics.get_topic("pkpd-export").body
    for col in ["time_h", "conc_mg_per_L", "gut_mg", "central_mg", "peripheral_mg", "effect_site_mg_per_L", "effect"]:
        in_code = any(col in Path(p).read_text(encoding="utf-8") for p in [PK, wf.__file__] if Path(p).exists())
        check(f"export column {col} exists in code", in_code)
        check(f"manual mentions {col}", col in export_body)
except Exception as exc:  # pragma: no cover
    check("drift import", False, repr(exc))

print(f"\n{passed} passed, {failed} failed")
sys.exit(1 if failed else 0)
