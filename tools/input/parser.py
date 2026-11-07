#!/usr/bin/env python3
"""Parse a YAML/JSON configuration or a Gendata.in file into the input model.

InputParser takes a schema - generated from the source tree by gendata_schema,
or loaded from a file with gendata_schema.load - and turns either format into
the model that InputGendata writes: a `Model` of key entries and block entries
in file order, where a value is a tagged tuple ('logical'|'integer'|'real'|
'string'|'array', data) and a block row a pair (tokens, comment). The 'pre' and
'comment' fields carry the comments of the file so that they survive the round
trip.

An unknown key is an error when parsing YAML (unless allow_unknown) and a
warning when unpacking a .in file, so that nothing is silently dropped. Former
names resolve to the present name of the key. A problem in the file raises
InputFormat; everything softer lands in the list attribute `warnings`.
"""

from __future__ import annotations

import difflib
import re
from dataclasses import dataclass, field
from typing import Any, Literal, TypedDict, cast

from gendata_schema import BlockInfo, KeyInfo, Schema, lookup, normalize

LOGICAL_TRUE = ('.true.', 'true', 't', 'yes', 'y')
LOGICAL_FALSE = ('.false.', 'false', 'f', 'no', 'n')

ARITHMETIC_RE = re.compile(r'[-+*/().0-9eEdD]+')  # what the solver's parser may evaluate

Scalar = bool | int | float | str
Token = tuple[Literal['logical', 'integer', 'real', 'string'], Scalar]
ArrayValue = tuple[Literal['array'], list[int | float]]
Value = Token | ArrayValue
BlockRow = tuple[list[Token], str | None]


class InputFormat(Exception):
    """A problem in the file being parsed: reported to the user without a traceback."""


class KeyEntry(TypedDict):
    """One key of the model: its value and the comments that came with it."""

    name: str  # present name, as the solver spells it
    entry: KeyInfo | None  # schema entry, None when the key is unknown to this version
    value: Value
    pre: list[str]  # comment lines above the key
    comment: str | None  # comment at the end of the key line


class BlockEntry(TypedDict):
    """One block of the model: rows of tokens instead of a value."""

    name: str
    entry: BlockInfo | None
    rows: list[BlockRow]
    pre: list[str]
    comment: str | None


Entry = KeyEntry | BlockEntry
OrderKind = Literal['keys', 'blocks']
OrderItem = tuple[OrderKind, str]


class Model(TypedDict):
    """A parsed input file: entries by normalized name, plus the order to write them in."""

    keys: dict[str, KeyEntry]
    blocks: dict[str, BlockEntry]
    order: list[OrderItem]
    trailing: list[str]  # comment lines after the last entry


def entry_of(model: Model, kind: OrderKind, key: str) -> Entry:
    """The entry of the item (kind, key) of model['order'], of either kind."""
    return model['keys'][key] if kind == 'keys' else model['blocks'][key]


def parse_token(tok: str) -> Token:
    """One word of a .in value into ('logical'|'integer'|'real'|'string', value).

    The solver's parser evaluates arithmetic in values, which is how 2/3 and the 0-0.01
    form of a negative number work; the same is done here for the words it accepts.
    """
    tl = tok.lower()
    if tl in LOGICAL_TRUE:
        return ('logical', True)
    if tl in LOGICAL_FALSE:
        return ('logical', False)
    if re.fullmatch(r'[-+]?\d+', tok):
        return ('integer', int(tok))
    try:
        return ('real', float(tok.replace('d', 'e').replace('D', 'E')))
    except ValueError:
        pass
    value = evaluate(tok)
    if value is not None:
        return ('integer', value) if isinstance(value, int) else ('real', float(value))
    return ('string', tok)


def evaluate(tok: str) -> int | float | None:
    """The value of an arithmetic expression of numbers and operators, or None.

    Only the restricted grammar of the solver's value parser is accepted, so a string
    value can never be executed; '2/3', '0-0.01' and '0-10.5/3.1+0.24' pass, 'Graphene'
    or 'AB-BA' do not.
    """
    if not ARITHMETIC_RE.fullmatch(tok):
        return None
    if not re.search(r'[-+*/]', tok.lstrip('+-')):
        return None
    try:
        value: int | float = eval(tok.replace('d', 'e').replace('D', 'E'), {'__builtins__': {}}, {})
    except Exception:
        return None
    return value


def split_comment(line: str) -> tuple[str, str | None]:
    """A .in line into (code, comment), as the solver reads it (InputProcess).

    A line whose first character is '!', '#' or '*' is a comment; elsewhere '!' or '#' starts
    a comment wherever it stands, with or without a space before it.
    """
    stripped = line.lstrip()
    if stripped[:1] in ('!', '#', '*'):
        return '', stripped[1:].strip().lstrip('!#*').strip()
    m = re.search(r'[!#]', line)
    if m:
        return line[: m.start()], line[m.start() + 1 :].strip().lstrip('!#').strip()
    return line, None


def solver_lines(text: str) -> list[tuple[str, str | None, int]]:
    """The logical lines of a .in text as (code, comment, number of the first physical line).

    As in the solver: a line whose code ends with '\\' continues on the next line (joined with
    a space); comment lines inside a continuation are skipped, and a blank line ends it.
    """
    out: list[tuple[str, str | None, int]] = []
    held: list[str] = []  # code of the lines of an open continuation
    held_comments: list[str] = []
    start = 0
    for number, line in enumerate(text.splitlines(), start=1):
        body, comment = split_comment(line)
        if held and not body.strip():
            if comment is not None:  # a comment line does not end a continuation
                held_comments.append(comment)
                continue
            out.append((' '.join(held), ' '.join(held_comments) or None, start))
            held, held_comments = [], []
        if not body.strip():
            out.append((body, comment, number))
            continue
        if not held:
            start = number
        if comment is not None:
            held_comments.append(comment)
        if body.rstrip().endswith('\\'):
            held.append(body.rstrip()[:-1].strip())
            continue
        held.append(body.strip())
        out.append((' '.join(held), ' '.join(held_comments) or None, start))
        held, held_comments = [], []
    if held:
        out.append((' '.join(held), ' '.join(held_comments) or None, start))
    return out


def clean_comment(text: str) -> str:
    """A comment token as plain text: no markers, no surrounding space."""
    return text.strip().lstrip('#').strip()


def comment_lines_of(text: str) -> list[str]:
    """The non-empty lines of a multi-line comment, each cleaned."""
    return [clean_comment(line) for line in text.splitlines() if line.strip()]


def comment_tokens(x: Any) -> list[str]:
    """The text of a ruamel comment slot: a token, a list of tokens, or nothing."""
    if x is None:
        return []
    return [t.value for t in (x if isinstance(x, list) else [x]) if t is not None]


def map_comments(m: Any) -> dict[Any, tuple[list[str], str | None]]:
    """{key: (pre_list, eol)} for a mapping read with ruamel in round-trip mode.

    ruamel stores a standalone comment in the slot of the *preceding* key (its value
    starts with a newline) and a comment before the first key in the map's own slot;
    both are carried to the key they precede.
    """
    out: dict[Any, tuple[list[str], str | None]] = {}
    carry: list[str] = []
    ca = getattr(m, 'ca', None)
    if ca is not None and getattr(ca, 'comment', None):
        for value in comment_tokens(ca.comment[1] if len(ca.comment) > 1 else None):
            carry += comment_lines_of(value)
    items = getattr(ca, 'items', {}) if ca else {}
    for key in m:
        slot = items.get(key)
        pre, eol = list(carry), None
        carry = []
        if slot:
            for value in comment_tokens(slot[0] if len(slot) > 0 else None):
                pre += comment_lines_of(value)
            for value in comment_tokens(slot[1] if len(slot) > 1 else None):
                pre += comment_lines_of(value)
            # slots 2 and 3 hold a token or (after a load of a file with several comment lines) a list of them
            for value in comment_tokens(slot[2] if len(slot) > 2 else None):
                if value.startswith('\n'):
                    carry += comment_lines_of(value)
                elif eol is None:
                    eol = clean_comment(value)
                else:
                    carry += comment_lines_of(value)
            # for a section (a mapping value) slot 3 holds the comment before its first key, which the
            # section's own map also carries: passing it on to the next key would write it twice
            if not isinstance(m[key], dict):
                for value in comment_tokens(slot[3] if len(slot) > 3 else None):
                    carry += comment_lines_of(value)
        out[key] = (pre, eol)
    return out


def seq_item_comment(seq: Any, idx: int) -> str | None:
    """The comment of one item of a sequence read or built in ruamel round-trip mode.

    For an item that is itself a sequence, ruamel keeps the comment written after the '-'
    in slot 0 (in memory) or slot 1 (after a load), as tokens.
    """
    try:
        slot = seq.ca.items.get(idx)
    except AttributeError:
        return None
    if not slot:
        return None
    parts = comment_tokens(slot[0] if len(slot) > 0 else None)
    parts += comment_tokens(slot[1] if len(slot) > 1 else None)
    if len(slot) > 2 and slot[2] is not None:
        parts.append(slot[2].value)
    return clean_comment(' '.join(p.strip() for p in parts)) if parts else None


def new_model() -> Model:
    """An empty model."""
    return {'keys': {}, 'blocks': {}, 'order': [], 'trailing': []}


@dataclass(slots=True)
class InputParser:
    """Parses a YAML/JSON configuration or a Gendata.in against `schema` into the model."""

    schema: Schema
    allow_unknown: bool = False
    warnings: list[str] = field(default_factory=list)
    # key -> (precedence of the value kept, name as written there); reset by each parse
    _key_rank: dict[str, tuple[int, str]] = field(default_factory=dict)

    def parse_yaml(self, data: Any) -> Model:
        """The nested {section: {key: value}} mapping (or flat full names) to the model."""
        model = new_model()
        self._key_rank = {}
        if not isinstance(data, dict):
            raise InputFormat('the input must be a mapping of sections and keys')
        for section, content in data.items():
            if not isinstance(content, dict):  # full key name at the top level
                pre, eol = map_comments(data).get(section, ([], None))
                self._add_entry(model, str(section), content, pre, eol)
                continue
            comments = map_comments(content)
            for key, value in content.items():
                pre, eol = comments.get(key, ([], None))
                self._add_entry(model, f'{section}.{key}', value, pre, eol)
        return model

    def _add_entry(self, model: Model, fullname: str, value: Any, pre: list[str], eol: str | None) -> None:
        kind, entry, canonical = lookup(self.schema, fullname)
        if kind is None:
            if not self.allow_unknown:
                candidates = list(self.schema['keys']) + list(self.schema['blocks'])
                near = difflib.get_close_matches(normalize(fullname), candidates, 1)
                hint = (
                    f' (did you mean {self.schema["keys"][near[0]]["name"]}?)'
                    if near and near[0] in self.schema['keys']
                    else ''
                )
                raise InputFormat(f'{fullname}: no key of this name in this version of the solver{hint}')
            self.warnings.append(
                f'{fullname}: no key of this name in this version of the solver, writing it unchanged'
            )
        elif normalize(canonical or fullname) != normalize(fullname):
            self.warnings.append(f'{fullname}: former name of {canonical}, writing {canonical}')
        name = canonical or fullname
        if kind == 'block':
            model['blocks'][normalize(name)] = {
                'name': name,
                'entry': cast('BlockInfo', entry),
                'rows': self._block_rows(value, fullname),
                'pre': pre,
                'comment': eol,
            }
            model['order'].append(('blocks', normalize(name)))
            return
        key_info = cast('KeyInfo | None', entry)
        action = self._repeat(name, fullname, key_info)
        if action == 'skip':
            return
        if action == 'replace':
            kept = model['keys'][normalize(name)]
            kept['value'] = self._scalar_value(value, key_info, fullname)
            kept['comment'] = eol
            return
        model['keys'][normalize(name)] = {
            'name': name,
            'entry': key_info,
            'value': self._scalar_value(value, key_info, fullname),
            'pre': pre,
            'comment': eol,
        }
        model['order'].append(('keys', normalize(name)))

    def _repeat(self, found: str, written: str, entry: KeyInfo | None) -> str:
        """'new' for the first line of a key; for a repeat, 'replace' or 'skip' as the solver would."""
        rank = self._precedence(written, found, entry)
        used = self._key_rank.get(normalize(found))
        if used is None:
            self._key_rank[normalize(found)] = (rank, written)
            return 'new'
        used_rank, used_as = used
        if rank < used_rank:
            self.warnings.append(
                f'{found}: given as {used_as} and as {written}; the value of {written} is used'
            )
            self._key_rank[normalize(found)] = (rank, written)
            return 'replace'
        if rank > used_rank:
            self.warnings.append(
                f'{found}: given as {used_as} and as {written}; the value of {used_as} is used'
            )
        else:
            self.warnings.append(f'{found}: appears twice, the first value is used')
        return 'skip'

    def _scalar_value(self, value: Any, entry: KeyInfo | None, fullname: str) -> Value:
        typ = entry['type'] if entry else 'unknown'
        if typ == 'logical':
            if isinstance(value, bool):
                return ('logical', value)
            if isinstance(value, str) and value.lower() in LOGICAL_TRUE + LOGICAL_FALSE:
                return ('logical', value.lower() in LOGICAL_TRUE)
            raise InputFormat(f'{fullname}: expected true or false, got {value!r}')
        if typ == 'integer':
            if isinstance(value, bool) or not isinstance(value, int):
                if isinstance(value, float) and value.is_integer():
                    self.warnings.append(
                        f'{fullname}: {value} is not an integer literal, writing {int(value)}'
                    )
                    return ('integer', int(value))
                number = evaluate(value) if isinstance(value, str) else None
                if number is not None:
                    return self._scalar_value(number, entry, fullname)
                raise InputFormat(f'{fullname}: expected an integer, got {value!r}')
            return ('integer', value)
        if typ == 'real':
            number = evaluate(value) if isinstance(value, str) else None
            if number is not None:  # '0-10.5/3.1'
                return ('real', float(number))
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise InputFormat(f'{fullname}: expected a number, got {value!r}')
            return ('real', float(value))
        if typ == 'string':
            if not isinstance(value, str):
                raise InputFormat(f'{fullname}: expected a string, got {value!r}')
            return ('string', value)
        if typ == 'array':
            return ('array', self._array_items(value, entry, fullname))
        if isinstance(value, bool):  # type not fixed by a literal default
            return ('logical', value)
        if isinstance(value, int):
            return ('integer', value)
        if isinstance(value, float):
            return ('real', value)
        if isinstance(value, list):
            return ('array', self._array_items(value, entry, fullname))
        if isinstance(value, str):
            return ('string', value)
        raise InputFormat(f'{fullname}: cannot write {value!r}')

    def _array_items(self, value: Any, entry: KeyInfo | None, fullname: str) -> list[int | float]:
        if not isinstance(value, list):
            raise InputFormat(f'{fullname}: expected a list of numbers')
        info = entry.get('array') if entry else None
        length = info['length'] if info else None
        if length and len(value) != length:
            raise InputFormat(f'{fullname}: expected {length} numbers, got {len(value)}')
        if any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in value):
            raise InputFormat(f'{fullname}: expected numbers, got {value!r}')
        return [float(x) if isinstance(x, float) else int(x) for x in value]

    def _block_rows(self, value: Any, fullname: str) -> list[BlockRow]:
        if not isinstance(value, list) or not value:
            raise InputFormat(f'{fullname}: a block needs a list of rows')
        if not isinstance(value[0], list):
            value = [value]  # a single row may be written flat
        rows: list[BlockRow] = []
        for idx, row in enumerate(value):
            if not isinstance(row, list):
                raise InputFormat(f'{fullname}: a block row must be a list of numbers')
            tokens: list[Token] = []
            for x in row:
                if isinstance(x, bool):
                    tokens.append(('logical', x))
                elif isinstance(x, str):
                    token = parse_token(x)
                    tokens.append(token if token[0] != 'string' else ('string', x))
                elif isinstance(x, int):
                    tokens.append(('integer', x))
                elif isinstance(x, float):
                    tokens.append(('real', x))
                else:
                    raise InputFormat(f'{fullname}: cannot write {x!r} in a block row')
            rows.append((tokens, seq_item_comment(value, idx)))
        return rows

    def parse_in(self, text: str, name: str = '<input>') -> Model:
        """The text of a Gendata.in into the model, with the solver's precedence for repeated keys.

        As in the solver (InputSearchLabel): the present (Section.Name) name of a key wins over its
        former names wherever it appears, then the former names in the order of the alias table;
        among lines of the same spelling, the first one wins.
        """
        model = new_model()
        lines = solver_lines(text)
        pending: list[str] = []  # comment lines before the next entry
        seen: set[str] = set()
        self._key_rank = {}
        i = 0
        while i < len(lines):
            body, comment, number = lines[i]
            words = body.split()
            if not words:
                if comment:
                    pending.append(comment)
                i += 1
                continue
            if words[0].lower() in ('&begin', '&end'):
                if words[0].lower() == '&end':
                    raise InputFormat(f'{name}, line {number}: "&end" without "&begin"')
                block_name = words[1]
                declared = int(words[2]) if len(words) > 2 and re.fullmatch(r'\d+', words[2]) else None
                rows: list[BlockRow] = []
                j = i + 1
                while j < len(lines):
                    body2, comment2, _ = lines[j]
                    if body2.split()[:1] == ['&end']:
                        break
                    if body2.strip():
                        rows.append(([parse_token(t) for t in body2.split()], comment2))
                    j += 1
                if j >= len(lines):
                    raise InputFormat(f'{name}, line {number}: block {block_name} is not closed with "&end"')
                if declared is not None and declared != len(rows):
                    self.warnings.append(
                        f'{block_name}: the block declares {declared} lines but contains {len(rows)}'
                    )
                kind, entry, canonical = lookup(self.schema, block_name)
                block = cast('BlockInfo', entry) if kind == 'block' else None
                if block is None:
                    self.warnings.append(f'{block_name}: not a block of this version of the solver')
                    canonical = block_name
                found = canonical or block_name
                if normalize(found) in seen:
                    self.warnings.append(f'{found}: appears twice, the first one is used')
                    i = j + 1
                    continue
                seen.add(normalize(found))
                model['blocks'][normalize(found)] = {
                    'name': found,
                    'entry': block,
                    'rows': rows,
                    'pre': pending,
                    'comment': comment,
                }
                model['order'].append(('blocks', normalize(found)))
                pending = []
                i = j + 1
                continue
            key, rest = words[0], body[len(words[0]) :].strip()
            kind, entry, canonical = lookup(self.schema, key)
            if kind is None:
                near = difflib.get_close_matches(
                    normalize(key), list(self.schema['keys']) + list(self.schema['blocks']), 1
                )
                hint = (
                    f' (did you mean {self.schema["keys"][near[0]]["name"]}?)'
                    if near and near[0] in self.schema['keys']
                    else ''
                )
                self.warnings.append(f'{key}: no key of this name in this version of the solver{hint}')
                canonical, entry = key, None
            elif canonical is not None and normalize(canonical) != normalize(key):
                self.warnings.append(f'{key}: former name of {canonical}, reading it as {canonical}')
            found = canonical or key
            key_info = cast('KeyInfo | None', entry) if kind == 'key' else None
            action = self._repeat(found, key, key_info)
            if action == 'replace':
                kept = model['keys'][normalize(found)]
                kept['value'] = self._parse_value(rest, key_info)
                kept['comment'] = comment
            if action != 'new':
                i += 1
                continue
            model['keys'][normalize(found)] = {
                'name': found,
                'entry': key_info,
                'value': self._parse_value(rest, key_info),
                'pre': pending,
                'comment': comment,
            }
            model['order'].append(('keys', normalize(found)))
            pending = []
            i += 1
        if not model['order']:
            raise InputFormat(f'{name}: no keys found')
        model['trailing'] = pending
        return model

    @staticmethod
    def _precedence(written: str, found: str, entry: KeyInfo | None) -> int:
        """0 for the present name of a key, 1, 2, ... for its former names in alias-table order."""
        if normalize(written) == normalize(found) or entry is None:
            return 0
        former = [normalize(x) for x in entry.get('formerNames', [])]
        return 1 + former.index(normalize(written)) if normalize(written) in former else 1 + len(former)

    @staticmethod
    def _parse_value(rest: str, entry: KeyInfo | None) -> Value:
        """The value part of a key line: the whole line for string keys, else one token per word."""
        if not rest:
            raise InputFormat('a key without a value')
        if entry and entry['type'] == 'string':
            return ('string', rest)
        tokens = [parse_token(t) for t in rest.split()]
        if len(tokens) == 1:
            return tokens[0]
        if all(t[0] in ('integer', 'real') for t in tokens):
            numbers = [t[1] for t in tokens]
            return ('array', [float(x) if isinstance(x, float) else int(x) for x in numbers])
        return ('string', rest)
