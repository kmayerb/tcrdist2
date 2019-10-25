import pytest
import pandas as pd
from collections import OrderedDict
from tcrdist import mappers

def test_generic_pandas_mapper():
    """
    Test that generic pandas mapping selects and renames columns
    """
    df = pd.DataFrame({"a":[1,2],"b":[3,4]})
    mapper = OrderedDict([('a', 'ardvark')])
    r = mappers.generic_pandas_mapper(df, mapper)
    expected = pd.DataFrame({"ardvark":[1,2]})
    assert r.equals(expected)

def test_generic_pandas_mapper_raises_assertion_error():
    """
    Test that a KeyError is raised if the mapper contains keys not in df.columns
    """
    df = pd.DataFrame({"a":[1,2],"b":[3,4]})
    mapper = OrderedDict([('c', 'canary')])
    with pytest.raises(KeyError):
        r = mappers.generic_pandas_mapper(df, mapper)
        expected = pd.DataFrame({"ardvark":[1,2]})
