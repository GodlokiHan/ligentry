"""Turn topics into pages for the Help window, and handle links and history (no Qt).

Links inside a page:

    topic:<id>              another manual page
    websearch:<source>:<q>  a search at a web source (the "look it up" links)
    https://...             a web page
"""
from __future__ import annotations

import html
from dataclasses import dataclass
from urllib.parse import quote, unquote

from ligentry.help.search import Hit
from ligentry.help.topics import SECTIONS, TOPICS, get_topic, has_topic, topics_in_section
from ligentry.help.web import QUICK_SOURCES, build_search_url, get_source, is_safe_url

HOME_ID = "home"
SEARCH_PREFIX = "search:"

_TABLE_ATTRS = 'border="1" cellspacing="0" cellpadding="4"'


def topic_link(topic_id: str, text: str) -> str:
    return f'<a href="topic:{topic_id}">{html.escape(text)}</a>'


def web_search_link(source_key: str, query: str, text: str | None = None) -> str:
    label = text if text is not None else get_source(source_key).name
    return f'<a href="websearch:{source_key}:{quote(query, safe="")}">{html.escape(label)}</a>'


def render_page(topic_id: str) -> str:
    """The full HTML page of one topic, with See also and Look it up sections."""
    t = get_topic(topic_id)
    parts = [f"<h2>{html.escape(t.title)}</h2>",
             f"<p><i>{html.escape(t.section)}</i></p>",
             t.body]
    if t.related:
        links = " &middot; ".join(topic_link(r, get_topic(r).title) for r in t.related)
        parts.append(f"<h4>See also</h4><p>{links}</p>")
    if t.web:
        rows = []
        for q in t.web:
            links = ", ".join(web_search_link(k, q) for k in QUICK_SOURCES)
            rows.append(f"<li><b>{html.escape(q)}</b>: {links}</li>")
        parts.append("<h4>Look it up on the web</h4><ul>" + "".join(rows) + "</ul>")
    return "\n".join(parts)


def render_home() -> str:
    """The contents page: every section with its topics."""
    parts = ["<h2>ligentry manual</h2>",
             "<p>Pick a page below, or type in the search box above. "
             "Press <kbd>F1</kbd> in any window to come back here at the right page.</p>"]
    for section in SECTIONS:
        items = "".join(
            f"<li>{topic_link(t.id, t.title)}: {html.escape(t.summary)}</li>"
            for t in topics_in_section(section)
        )
        parts.append(f"<h3>{html.escape(section)}</h3><ul>{items}</ul>")
    return "\n".join(parts)


def render_search_results(query: str, hits: list[Hit]) -> str:
    q = html.escape(" ".join(query.split()))
    if not hits:
        web = ", ".join(web_search_link(k, query) for k in QUICK_SOURCES) if query.strip() else ""
        extra = f"<p>Search the web instead: {web}</p>" if web else ""
        return (f"<h2>Search results</h2><p>Nothing in the manual matches <b>{q}</b>. "
                f"Try fewer or simpler words, or open the {topic_link(HOME_ID, 'contents')}.</p>{extra}")
    items = "".join(
        f"<li><b>{topic_link(h.topic.id, h.topic.title)}</b> <i>({html.escape(h.topic.section)})</i><br>"
        f"{html.escape(h.snippet)}</li>"
        for h in hits
    )
    web = ", ".join(web_search_link(k, query) for k in QUICK_SOURCES)
    return (f"<h2>Search results</h2><p>{len(hits)} page(s) match <b>{q}</b>:</p><ul>{items}</ul>"
            f"<h4>Not found what you need?</h4><p>Search the web for it: {web}</p>")


def render(page_id: str, hits: list[Hit] | None = None) -> str:
    """Page for a history entry: 'home', 'search:<query>' (needs hits) or a topic id."""
    if page_id == HOME_ID:
        return render_home()
    if page_id.startswith(SEARCH_PREFIX):
        return render_search_results(page_id[len(SEARCH_PREFIX):], hits or [])
    return render_page(page_id)


@dataclass(frozen=True)
class Link:
    kind: str  # "topic", "web", "blocked" or "unknown"
    value: str  # topic id, address, or the original text
    label: str = ""


def classify_link(href: str) -> Link:
    """Decide what clicking a link in a manual page should do."""
    href = href.strip()
    if href.startswith("topic:"):
        tid = href[len("topic:"):]
        if tid == HOME_ID or has_topic(tid):
            return Link("topic", tid)
        return Link("unknown", href)
    if href.startswith("websearch:"):
        try:
            _, key, q = href.split(":", 2)
            url = build_search_url(key, unquote(q))
        except ValueError:
            return Link("unknown", href)
        return Link("web", url, f"{get_source(key).name}: {unquote(q)}")
    if href.lower().startswith(("http://", "https://")):
        return Link("web", href, href) if is_safe_url(href) else Link("blocked", href)
    return Link("unknown", href)


class History:
    """Back and forward through pages, like a browser (QTextBrowser.setHtml keeps none)."""

    def __init__(self, start: str = HOME_ID):
        self._items = [start]
        self._pos = 0

    @property
    def current(self) -> str:
        return self._items[self._pos]

    @property
    def can_back(self) -> bool:
        return self._pos > 0

    @property
    def can_forward(self) -> bool:
        return self._pos < len(self._items) - 1

    def visit(self, page_id: str) -> str:
        """Go to a new page; forward history is dropped. Revisiting the current page changes nothing."""
        if page_id != self.current:
            del self._items[self._pos + 1:]
            self._items.append(page_id)
            self._pos += 1
        return self.current

    def back(self) -> str:
        if self.can_back:
            self._pos -= 1
        return self.current

    def forward(self) -> str:
        if self.can_forward:
            self._pos += 1
        return self.current


def all_topic_ids() -> list[str]:
    return [t.id for t in TOPICS]
