#!/usr/bin/env python3
"""Convert between a YAML/JSON input file and the solver's Gendata.in format.

The YAML file mirrors the Section.Name form of the input keys:

    run:
      description: pristine graphene, bands
      typeOfSystem: Graphene
    structure:
      cellSize: 6
      superCell: 1
    tb:
      hopping: 3.1
    moire:
      twistAngle: 1.05
      potential: true
    bands:
      numPoints: 40
      path:                       # the Bands.Path block
        - [0.0, 0.0, 0.0]         # Gamma
        - [0.5, 0.0, 0.0]         # M

InputParser checks the keys and types against the schema of this checkout
(gendata_schema), so a misspelt key is an error here instead of a silent
default in the solver; InputGendata writes the model back as Gendata.in, YAML
or JSON. The conversion also fixes the two formats' quirks: booleans become
.true., negative numbers become the 0-0.01 form, and comments survive the
round trip.

    python3 tools/input/gendata.py pack config.yaml -o Gendata.in
    python3 tools/input/gendata.py unpack Gendata.in -o config.yaml
    python3 tools/input/gendata.py dump-schema schema.json

With --schema FILE the conversion uses a stored schema instead of generating
one from the source tree, which is how an input is checked against another
version of the code.
"""

from __future__ import annotations

import functools
import io
import json
import os
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any, cast

import click
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gendata_schema import Schema, generate, load, save  # noqa: E402
from parser import (  # noqa: E402
    InputFormat,
    InputParser,
    Model,
    OrderItem,
    OrderKind,
    Token,
    Value,
    entry_of,
)

SECTIONS = [
    'Run',
    'Structure',
    'TB',
    'Moire',
    'moireCDW',
    'moireBLShift',
    'GBN',
    'BNBN',
    'Calculate',
    'Diag',
    'TAPW',
    'Bands',
    'Kubo',
    'Spectral',
    'MagField',
    'Haldane',
    'HaldPhase',
    'SOC',
    'Interlayer',
    'Intralayer',
    'Hopping',
    'Potential',
    'Disorder',
    'Neigh',
    'Output',
    'Strain',
    'Hubbard',
    'Hybrid',
    'Interface',
    'Random',
    'Input',
]

KEY_WIDTH = 24


def fmt_token(token: Token) -> str:
    """A parsed token in the .in spelling: .true., and the 0-5 / 0-0.01 form of negative numbers.

    The solver reads integers with the same expression parser as reals, so a bare -5 stops a
    run with "Bad real number" just like -0.01.
    """
    kind, value = token
    if kind == 'logical':
        return '.true.' if value else '.false.'
    if kind == 'integer':
        s = str(value)
        return '0' + s if s.startswith('-') else s
    if kind == 'real':
        s = repr(float(value))
        return '0' + s if s.startswith('-') else s  # the 0-0.01 form
    return str(value)


def fmt_value(value: Value) -> str:
    """A parsed value in the .in spelling: one word, or words side by side for an array."""
    if value[0] == 'array':
        items = cast('list[int | float]', value[1])
        return '  '.join(fmt_token(('real' if isinstance(x, float) else 'integer', x)) for x in items)
    return fmt_token(cast('Token', value))


def section_of(name: str) -> str:
    """The section of a key name: 'Bands.Path' gives 'Bands', 'MoireTwistAngle' gives ''."""
    return name.split('.')[0] if '.' in name else ''


@dataclass(slots=True)
class InputGendata:
    """Writes a parsed model (see parser.InputParser) out as Gendata.in, YAML or JSON."""

    schema: Schema

    def to_infile(self, model: Model) -> str:
        """The model to .in text: keys grouped by section, blocks after the scalars."""
        groups: dict[str, list[OrderItem]] = {}
        for kind, key in model['order']:
            groups.setdefault(section_of(entry_of(model, kind, key)['name']), []).append((kind, key))

        def rank(section: str) -> tuple[int, str]:
            return (SECTIONS.index(section) if section in SECTIONS else len(SECTIONS), section.lower())

        out: list[str] = []
        for section in sorted(groups, key=rank):
            body: list[str] = []
            for kind, key in sorted(
                groups[section],
                key=lambda ok: (ok[0] == 'blocks', entry_of(model, ok[0], ok[1])['name'].lower()),
            ):
                entry = entry_of(model, kind, key)
                for pre in entry['pre']:
                    body.append(f'# {pre}')
                if kind == 'keys':
                    value = model['keys'][key]['value']
                    line = f'{entry["name"]} {fmt_value(value)}'
                    if len(line) < KEY_WIDTH:
                        line = f'{entry["name"]:<{KEY_WIDTH}}{fmt_value(value)}'
                    body.append(line + (f'   # {entry["comment"]}' if entry['comment'] else ''))
                else:
                    rows = model['blocks'][key]['rows']
                    body.append(
                        f'&begin {entry["name"]} {len(rows)}'
                        + (f'   # {entry["comment"]}' if entry['comment'] else '')
                    )
                    for row, rowcomment in rows:
                        body.append(
                            '  '.join(fmt_token(t) for t in row)
                            + (f'   # {rowcomment}' if rowcomment else '')
                        )
                    body.append(f'&end {entry["name"]}')
            if out:
                out.append('')
            out += body
        for trailing in model['trailing']:
            out.append(f'# {trailing}')
        return '\n'.join(out) + '\n'

    def to_yaml(self, model: Model) -> CommentedMap:
        """The model to nested CommentedMaps; comments of the model are carried over."""

        def sortkey(ok: OrderItem) -> tuple[int, str, bool, str]:
            name = entry_of(model, ok[0], ok[1])['name']
            section = section_of(name)
            return (
                SECTIONS.index(section) if section in SECTIONS else len(SECTIONS),
                section.lower(),
                ok[0] == 'blocks',
                name.lower(),
            )

        doc = CommentedMap()
        for kind, key in sorted(model['order'], key=sortkey):
            entry = entry_of(model, kind, key)
            name = entry['name']
            if '.' in name:
                section, rest = name.split('.', 1)
                target: CommentedMap = doc.setdefault(section, CommentedMap())
            else:
                target, rest = doc, name
            ykey = rest[0].lower() + rest[1:] if rest else rest
            target[ykey] = self._plain_value(model, kind, key)
            if entry['comment']:
                target.yaml_add_eol_comment(entry['comment'], ykey)
            if entry['pre']:
                target.yaml_set_comment_before_after_key(ykey, before='\n'.join(entry['pre']))
        return doc

    def to_json(self, model: Model) -> dict[str, Any]:
        """The model to a plain {key: value} dict, block names included; no comments."""
        data: dict[str, Any] = {}
        for kind, key in model['order']:
            name = entry_of(model, kind, key)['name']
            if kind == 'blocks':
                data[name] = [[t[1] for t in row] for row, _ in model['blocks'][key]['rows']]
            else:
                data[name] = self._plain_value(model, kind, key)
        return data

    @staticmethod
    def _plain_value(model: Model, kind: OrderKind, key: str) -> Any:
        """The value of an entry as a plain YAML/JSON object: a list of rows for a block."""
        if kind == 'blocks':
            seq = CommentedSeq()
            for row, comment in model['blocks'][key]['rows']:
                seq.append([t[1] for t in row])
                if comment:
                    seq.yaml_add_eol_comment(comment, len(seq) - 1)
            return seq
        value = model['keys'][key]['value']
        if value[0] == 'array':
            return list(cast('list[int | float]', value[1]))
        return value[1]


def yaml_writer() -> YAML:
    """A round-trip YAML writer: comments survive, and long lines are not folded."""
    y = YAML()
    y.width = 4096
    y.preserve_quotes = True
    return y


def get_schema(schema_file: str | None) -> Schema:
    """The stored schema, or the one generated from the source tree of this checkout."""
    return load(schema_file) if schema_file else generate()


def catch_input_format(f: Callable[..., None]) -> Callable[..., None]:
    """Report an InputFormat raised by a command as a click error: message, exit 1, no traceback."""

    @functools.wraps(f)
    def wrapper(*args: Any, **kwargs: Any) -> None:
        try:
            f(*args, **kwargs)
        except InputFormat as e:
            raise click.ClickException(str(e)) from None

    return wrapper


def report(messages: Sequence[str], prefix: str = 'warning') -> None:
    """Print the parser's warnings to stderr, one per line."""
    for message in messages:
        click.echo(f'{prefix}: {message}', err=True)


@click.group()
def main() -> None:
    """Convert between YAML/JSON and the solver's Gendata.in format."""


@main.command()
@click.argument('config', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(), help='write to this file (default: stdout)')
@click.option(
    '--schema',
    'schema_file',
    type=click.Path(exists=True),
    help='stored schema to validate against (default: generated from the source tree)',
)
@click.option('--allow-unknown', is_flag=True, help='write keys unknown to this version unchanged')
@click.option('--check', is_flag=True, help='only validate, write nothing')
@catch_input_format
def pack(config: str, output: str | None, schema_file: str | None, allow_unknown: bool, check: bool) -> None:
    """Convert a YAML/JSON input to Gendata.in."""
    schema = get_schema(schema_file)
    ext = os.path.splitext(config)[1].lower()
    if ext == '.json':
        with open(config) as f:
            data: Any = json.load(f)
    elif ext in ('.yaml', '.yml'):
        with open(config) as f:
            data = yaml_writer().load(f)
    else:
        raise InputFormat(f'{config}: not a .yaml, .yml or .json file')
    parser = InputParser(schema, allow_unknown=allow_unknown)
    model = parser.parse_yaml(data)
    text = InputGendata(schema).to_infile(model)
    if check:
        click.echo(
            f'{config}: {len(model["keys"])} keys and {len(model["blocks"])} blocks, '
            f'{len(parser.warnings)} warning(s)'
        )
        return
    report(parser.warnings)
    if output:
        with open(output, 'w') as f:
            f.write(text)
        click.echo(
            f'wrote {output}: {len(model["keys"])} keys and {len(model["blocks"])} blocks '
            f'({len(parser.warnings)} warning(s))'
        )
    else:
        click.echo(text, nl=False)


@main.command()
@click.argument('infile', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(), help='write to this file (default: stdout)')
@click.option(
    '--schema',
    'schema_file',
    type=click.Path(exists=True),
    help='stored schema (default: generated from the source tree)',
)
@click.option('--format', 'fmt', type=click.Choice(['yaml', 'json']), default='yaml')
@catch_input_format
def unpack(infile: str, output: str | None, schema_file: str | None, fmt: str) -> None:
    """Convert a Gendata.in to YAML (or JSON with --format json)."""
    schema = get_schema(schema_file)
    parser = InputParser(schema)
    with open(infile) as f:
        model = parser.parse_in(f.read(), name=infile)
    gendata = InputGendata(schema)
    if fmt == 'json':
        text = json.dumps(gendata.to_json(model), indent=1) + '\n'
    else:
        stream = io.StringIO()
        yaml_writer().dump(gendata.to_yaml(model), stream)
        text = stream.getvalue()
    report(parser.warnings)
    if output:
        with open(output, 'w') as f:
            f.write(text)
        click.echo(
            f'wrote {output}: {len(model["keys"])} keys and {len(model["blocks"])} blocks '
            f'({len(parser.warnings)} warning(s))'
        )
    else:
        click.echo(text, nl=False)


@main.command('dump-schema')
@click.argument('output', type=click.Path(), required=False)
def dump_schema(output: str | None) -> None:
    """Write the input schema of this source tree."""
    schema = generate()
    if output:
        save(schema, output)
        click.echo(
            f'wrote {output}: {len(schema["keys"])} keys, {len(schema["blocks"])} blocks, '
            f'grabnes {schema["grabnes"]["version"]} ({schema["grabnes"]["commit"]})'
        )
    else:
        click.echo(json.dumps(schema, indent=1, sort_keys=True))


if __name__ == '__main__':
    main()
