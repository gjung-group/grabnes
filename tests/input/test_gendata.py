"""Tests for the YAML/JSON to Gendata.in converter (tools/input/parser.py, gendata.py)."""

from __future__ import annotations

import glob
import os
import re
from pathlib import Path
from typing import Any, cast

import pytest
from click.testing import CliRunner
from gendata import InputGendata, fmt_token, main
from gendata_schema import KeyInfo, Lookup, Schema, lookup
from parser import InputFormat, InputParser, Model, parse_token
from ruamel.yaml import YAML

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
EXAMPLE_INPUTS = sorted(glob.glob(os.path.join(ROOT, 'examples', '**', 'Gendata*.in'), recursive=True))


def comparable(model: Model) -> tuple[dict[str, Any], dict[str, Any]]:
    """The model reduced to values only, with floats rounded: what a round trip must preserve."""

    def number(x: object) -> object:
        return round(x, 9) if isinstance(x, float) else x

    keys: dict[str, Any] = {}
    for name, entry in model['keys'].items():
        kind, data = entry['value']
        if isinstance(data, list):
            keys[name] = (kind, [number(x) for x in data])
        else:
            keys[name] = (kind, number(data))
    blocks = {
        name: [[(kind, number(v)) for kind, v in row] for row, _ in entry['rows']]
        for name, entry in model['blocks'].items()
    }
    return keys, blocks


def roundtrip_in_in(schema: Schema, text: str) -> tuple[Model, Model]:
    """Parse .in text, write it out and parse it again."""
    parser = InputParser(schema, allow_unknown=True)
    model = parser.parse_in(text, name='test')
    text2 = InputGendata(schema).to_infile(model)
    return model, InputParser(schema, allow_unknown=True).parse_in(text2, name='again')


def key_entry(result: Lookup) -> KeyInfo:
    """The key entry of a lookup result; a block or an unknown name fails the test."""
    kind, entry, _ = result
    assert kind == 'key' and entry is not None
    return cast('KeyInfo', entry)


class TestTokens:
    def test_logical(self) -> None:
        for word in ('.true.', 'true', 'True', 'yes', 't'):
            assert parse_token(word) == ('logical', True)
        assert fmt_token(('logical', True)) == '.true.'
        assert fmt_token(('logical', False)) == '.false.'

    def test_negative_zero_prefix(self) -> None:
        assert parse_token('0-0.01') == ('real', -0.01)
        assert fmt_token(parse_token('0-0.01')) == '0-0.01'

    def test_bare_negative_is_accepted(self) -> None:
        assert parse_token('-0.01') == ('real', -0.01)
        assert fmt_token(parse_token('-0.01')) == '0-0.01'

    def test_fraction_and_d_exponent(self) -> None:
        assert parse_token('2/3') == ('real', pytest.approx(2 / 3))
        assert parse_token('135d0') == ('real', 135.0)

    def test_string(self) -> None:
        assert parse_token('Graphene') == ('string', 'Graphene')


class TestSchema:
    def test_keys_present(self, schema: Schema) -> None:
        assert len(schema['keys']) > 700
        assert key_entry(lookup(schema, 'Moire.TwistAngle'))['type'] == 'real'

    def test_former_name_resolves(self, schema: Schema) -> None:
        kind, _, canonical = lookup(schema, 'MoireTwistAngle')
        assert (kind, canonical) == ('key', 'Moire.TwistAngle')

    def test_case_and_underscore_insensitive(self, schema: Schema) -> None:
        assert lookup(schema, 'moire_TWIST-angle').name == 'Moire.TwistAngle'

    def test_array_key(self, schema: Schema) -> None:
        entry = key_entry(lookup(schema, 'Diag.KGrid'))
        assert entry['type'] == 'array'
        array = entry.get('array')
        assert array is not None and array['length'] == 3

    def test_block_and_pattern(self, schema: Schema) -> None:
        assert lookup(schema, 'Bands.Path').kind == 'block'
        assert lookup(schema, 'TB.Hoppings.2').kind == 'block'


class TestPack:
    def test_types_and_format(self, schema: Schema) -> None:
        model = InputParser(schema).parse_yaml(
            {
                'run': {'typeOfSystem': 'Graphene', 'description': 'a run'},
                'structure': {'cellSize': 6},
                'tb': {'hopping': 3.1},
                'moire': {'twistAngle': -1.05, 'twisted': True},
                'diag': {'kGrid': [4, 4, 1]},
                'bands': {'path': [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]},
            }
        )
        text = InputGendata(schema).to_infile(model)
        assert re.search(r'^Run\.TypeOfSystem\s+Graphene$', text, re.M)
        assert re.search(r'^Structure\.CellSize\s+6$', text, re.M)
        assert re.search(r'^Moire\.TwistAngle\s+0-1\.05$', text, re.M)
        assert re.search(r'^Moire\.Twisted\s+\.true\.$', text, re.M)
        assert re.search(r'^Diag\.KGrid\s+4\s+4\s+1$', text, re.M)
        assert '&begin Bands.Path 2' in text and '&end Bands.Path' in text

    def test_unknown_key_raises(self, schema: Schema) -> None:
        with pytest.raises(InputFormat, match='no key of this name'):
            InputParser(schema).parse_yaml({'moire': {'twistAngl': 1.05}})

    def test_close_match_suggested(self, schema: Schema) -> None:
        with pytest.raises(InputFormat, match='Moire.TwistAngle'):
            InputParser(schema).parse_yaml({'moire': {'twistAngleX': 1.0}})

    def test_wrong_type_raises(self, schema: Schema) -> None:
        with pytest.raises(InputFormat, match='expected an integer'):
            InputParser(schema).parse_yaml({'structure': {'cellSize': 'six'}})
        with pytest.raises(InputFormat, match='expected true or false'):
            InputParser(schema).parse_yaml({'moire': {'twisted': 1}})

    def test_array_length_checked(self, schema: Schema) -> None:
        with pytest.raises(InputFormat, match='expected 3 numbers'):
            InputParser(schema).parse_yaml({'diag': {'kGrid': [4, 4]}})

    def test_former_name_renamed_with_warning(self, schema: Schema) -> None:
        parser = InputParser(schema)
        model = parser.parse_yaml({'moireTwistAngle': 1.05})  # old section-less name, top level
        text = InputGendata(schema).to_infile(model)
        assert 'Moire.TwistAngle' in text and 'MoireTwistAngle' not in text
        assert any('former name' in w for w in parser.warnings)

    def test_unknown_key_allowed_when_asked(self, schema: Schema) -> None:
        parser = InputParser(schema, allow_unknown=True)
        model = parser.parse_yaml({'legacy': {'someOldKey': True}})
        assert re.search(r'^legacy\.someOldKey\s+\.true\.$', InputGendata(schema).to_infile(model), re.M)
        assert any('no key of this name' in w for w in parser.warnings)

    def test_comments_survive(self, schema: Schema) -> None:
        doc = YAML().load('moire:\n  # a comment above\n  twisted: true   # and beside\n')
        model = InputParser(schema).parse_yaml(doc)
        text = InputGendata(schema).to_infile(model)
        assert '# a comment above' in text
        assert 'Moire.Twisted           .true.   # and beside' in text


class TestParseIn:
    def test_first_value_wins(self, schema: Schema) -> None:
        parser = InputParser(schema)
        model = parser.parse_in('TB.Hopping 3.1\nTB.Hopping 2.9\n', name='test')
        assert model['keys']['tb.hopping']['value'] == ('real', 3.1)
        assert any('appears twice' in w for w in parser.warnings)

    def test_comments_as_in_the_solver(self, schema: Schema) -> None:
        # '!' and '#' start a comment anywhere; '!', '#' and '*' in the first column mark a comment line
        text = '* star line\n! bang line\nTB.Hopping 3.1!x\nKubo.Epsilon 0.002#y\n'
        parser = InputParser(schema)
        model = parser.parse_in(text, name='test')
        assert model['keys']['tb.hopping']['value'] == ('real', 3.1)
        assert model['keys']['kubo.epsilon']['value'] == ('real', 0.002)
        assert set(model['keys']) == {'tb.hopping', 'kubo.epsilon'}
        assert not parser.warnings

    def test_continuation_line(self, schema: Schema) -> None:
        # a trailing backslash joins the next line, also between key and value (the solver ran 60 steps)
        text = 'Kubo.RecursionNumber \\\n   60   # sixty\nKubo.Epsilon \\\n# inside\n 0.002\n'
        parser = InputParser(schema)
        model = parser.parse_in(text, name='test')
        assert model['keys']['kubo.recursionnumber']['value'] == ('integer', 60)
        assert model['keys']['kubo.recursionnumber']['comment'] == 'sixty'
        assert model['keys']['kubo.epsilon']['value'] == ('real', 0.002)
        assert not parser.warnings

    def test_present_name_wins_over_former_name(self, schema: Schema) -> None:
        # the solver uses Kubo.RecursionNumber wherever it is (checked: 60 recursion steps for both orders)
        for text in (
            'RecursionNumber 40\nKubo.RecursionNumber 60\n',
            'Kubo.RecursionNumber 60\nRecursionNumber 40\n',
        ):
            parser = InputParser(schema)
            model = parser.parse_in(text, name='test')
            assert model['keys']['kubo.recursionnumber']['value'] == ('integer', 60)
            assert any('the value of Kubo.RecursionNumber is used' in w for w in parser.warnings)

    def test_present_name_wins_in_yaml(self, schema: Schema) -> None:
        for data in (
            {'RecursionNumber': 40, 'kubo': {'recursionNumber': 60}},
            {'kubo': {'recursionNumber': 60}, 'RecursionNumber': 40},
        ):
            model = InputParser(schema).parse_yaml(data)
            assert model['keys']['kubo.recursionnumber']['value'] == ('integer', 60)
            assert model['order'].count(('keys', 'kubo.recursionnumber')) == 1

    def test_unclosed_block_raises(self, schema: Schema) -> None:
        with pytest.raises(InputFormat, match='not closed'):
            InputParser(schema).parse_in('&begin Bands.Path 2\n0.0 0.0 0.0\n', name='test')

    def test_end_without_begin_raises(self, schema: Schema) -> None:
        with pytest.raises(InputFormat, match='without "&begin"'):
            InputParser(schema).parse_in('&end Bands.Path\n', name='test')

    def test_comments_kept(self, schema: Schema) -> None:
        model = InputParser(schema).parse_in('# above\nTB.Hopping 3.1   # beside\n', name='test')
        entry = model['keys']['tb.hopping']
        assert entry['pre'] == ['above'] and entry['comment'] == 'beside'

    def test_string_value_is_the_rest_of_the_line(self, schema: Schema) -> None:
        model = InputParser(schema).parse_in('Run.Description two words here\n', name='test')
        assert model['keys']['run.description']['value'] == ('string', 'two words here')


class TestRoundTrip:
    @pytest.mark.parametrize('path', EXAMPLE_INPUTS, ids=lambda p: os.path.relpath(p, ROOT))
    def test_example_in_in(self, schema: Schema, path: str) -> None:
        with open(path) as f:
            model, again = roundtrip_in_in(schema, f.read())
        assert comparable(model) == comparable(again)

    @pytest.mark.parametrize('path', EXAMPLE_INPUTS, ids=lambda p: os.path.relpath(p, ROOT))
    def test_example_in_yaml_in(self, schema: Schema, path: str) -> None:
        parser = InputParser(schema)
        with open(path) as f:
            model = parser.parse_in(f.read(), name=path)
        yam = InputGendata(schema).to_yaml(model)
        parser2 = InputParser(schema, allow_unknown=True)
        assert comparable(model) == comparable(parser2.parse_yaml(yam))

    def test_in_yaml_in_keeps_comments(self, schema: Schema) -> None:
        text = '# a note\nTB.Hopping 3.1   # g0\n\n&begin Bands.Path 1\n0.0  0.0  0.0   # Gamma\n&end Bands.Path\n'
        model = InputParser(schema).parse_in(text, name='test')
        yam = InputGendata(schema).to_yaml(model)
        back = InputParser(schema).parse_yaml(yam)
        text2 = InputGendata(schema).to_infile(back)
        assert '# a note' in text2 and '# g0' in text2 and '# Gamma' in text2


class TestCli:
    def test_pack_unpack(self, schema_file: str, tmp_path: Path) -> None:
        config = tmp_path / 'config.yaml'
        config.write_text('run:\n  typeOfSystem: Graphene\ntb:\n  hopping: 3.1\ndiag:\n  calc: true\n')
        infile = tmp_path / 'Gendata.in'
        result = CliRunner().invoke(main, ['pack', str(config), '--schema', schema_file, '-o', str(infile)])
        assert result.exit_code == 0, result.output
        assert 'TB.Hopping' in infile.read_text()
        back = tmp_path / 'back.yaml'
        result = CliRunner().invoke(main, ['unpack', str(infile), '--schema', schema_file, '-o', str(back)])
        assert result.exit_code == 0, result.output
        assert 'typeOfSystem: Graphene' in back.read_text()

    def test_unknown_key_fails_with_suggestion(self, schema_file: str, tmp_path: Path) -> None:
        config = tmp_path / 'config.yaml'
        config.write_text('moire:\n  twistAngleX: 1.05\n')
        result = CliRunner().invoke(main, ['pack', str(config), '--schema', schema_file])
        assert result.exit_code == 1
        assert 'Moire.TwistAngle' in result.output

    def test_check_writes_nothing(self, schema_file: str, tmp_path: Path) -> None:
        config = tmp_path / 'config.yaml'
        config.write_text('tb:\n  hopping: 3.1\n')
        result = CliRunner().invoke(main, ['pack', '--check', str(config), '--schema', schema_file])
        assert result.exit_code == 0
        assert '1 keys' in result.output
        assert list(tmp_path.glob('Gendata.in')) == []

    def test_dump_schema(self, tmp_path: Path) -> None:
        out = tmp_path / 'schema.json'
        result = CliRunner().invoke(main, ['dump-schema', str(out)])
        assert result.exit_code == 0
        assert 'keys' in out.read_text()
