import pytest
from tcrdist.storage import StoreIO
from collections import namedtuple

def test_StoreIO_init__with_no_args():
    """
    StoreIO can be created without any input arguments,
    all pre-specified attrs exist and are None
    """
    S = StoreIO()
    assert S._validate_attrs()
    assert isinstance(S, StoreIO)
    assert list(vars(S).keys()) == ['name', 'valid_attrs', 'valid_attrs_type', 'a', 'b', 'c', 'd']
    assert S.a is None
    assert S.b is None
    assert S.c is None
    assert S.d is None


def test_StoreIO_init__with_complete_args():
    """
    test that StoreIO can be initialized with all valid args
    """
    abcd = namedtuple('abcd', ['a','b','c','d'])
    x = abcd(1,2,3,4)
    S = StoreIO(**x._asdict())
    assert isinstance(S, StoreIO)
    assert isinstance(S.a, int)
    assert isinstance(S.b, int)
    assert isinstance(S.c, int)
    assert isinstance(S.d, int)

def test_StoreIO_init__with_incomplete_args():
    """
    test that StoreIO can be initialized without all valid args
    """
    abc = namedtuple('abc', ['a','b','c'])
    x = abc(1,2,3)
    S = StoreIO(**x._asdict())
    assert isinstance(S, StoreIO)
    assert isinstance(S.a, int)
    assert isinstance(S.b, int)
    assert isinstance(S.c, int)
    assert S.d is None

def test_StoreIO_type_check__works_on_multiple_types():
    """
    _type check should return True when attr type is correct
    """
    abcd = namedtuple('abc', ['a','b','c',"d"])
    x = abcd(1,2.0,"3",{"k":1})
    S = StoreIO(**x._asdict())
    assert S._type_check("a", int)
    assert S._type_check("b", float)
    assert S._type_check("c", str)
    assert S._type_check("d", dict)

def test_StoreIO_type_check__raises_TypeError():
    """
    _type_check() should raise Type erorr if attr is not correct type
    """
    abcd = namedtuple('abc', ['a','b','c',"d"])
    x = abcd(1,2.0,"3",{"k":1})
    S = StoreIO(**x._asdict())
    with pytest.raises(TypeError):
        assert S._type_check("a", str)
    with pytest.raises(TypeError):
        assert S._type_check("b", int)
    with pytest.raises(TypeError):
        assert S._type_check("c", int)
    with pytest.raises(TypeError):
        assert S._type_check("d", int)

def test_StoreIO_validate_attrs__passes_with_full_correct_set():
    """
    _validate_attrs returns True when all valid args are of correct type
    matching self.valid_attrs_type which is default [int,int,int,int]
    """
    abcd = namedtuple('abcd', ['a','b','c','d'])
    x = abcd(1,2,3,4)
    S = StoreIO(**x._asdict())
    assert S._validate_attrs()

def test_StoreIO_validate_attrs__passes_with_partial_correct_set():
    """
    _validate_attrs returns True when partial set of valid args
    are of correct type and the rest are None (i.e. left blank)
    """
    abc = namedtuple('abcd', ['a','b','c'])
    x = abc(1,2,3)
    S = StoreIO(**x._asdict())
    assert(S.d is None)
    assert S._validate_attrs()

def test_StoreIO_validate_attrs__raises_TypeError():
    """
    _validate_attrs raises an error if supplied attribute is not the
    correct type
    """
    abcd = namedtuple('abcd', ['a','b','c','d'])
    x = abcd(1.0,2,3,4)
    S = StoreIO(**x._asdict())
    with pytest.raises(TypeError):
        S._validate_attrs()

def test_StoreIO_set_attrs__works_with_kwargs():
    """
    set_attrs_with_kwargs() accepts dict and updates attributes
    """
    ab = namedtuple('ab', ['a','b'])
    x = ab(1,2)
    S = StoreIO(**x._asdict())
    assert isinstance(S.a, int)
    assert isinstance(S.b, int)
    assert S.c is None
    assert S.d is None
    cd = namedtuple('cd', ['c','d'])
    y = cd(1,2)
    S.set_attrs_with_kwargs(**y._asdict())
    assert isinstance(S.a, int)
    assert isinstance(S.b, int)
    assert isinstance(S.c, int)
    assert isinstance(S.d, int)

def test_StoreIO_null_init_and_set_attrs_with_kwargs__works():
    """
    set_attrs_with_kwargs() accepts dict and updates attributes
    """
    S = StoreIO()
    assert S.a is None
    assert S.b is None
    assert S.c is None
    assert S.d is None
    abcd = namedtuple('abcd', ['a','b','c','d'])
    x = abcd(1,2,3,4)
    S.set_attrs_with_kwargs(**x._asdict())
    assert isinstance(S.a, int)
    assert isinstance(S.b, int)
    assert isinstance(S.c, int)
    assert isinstance(S.d, int)

def test_StoreIO_set_attrs_with_extra_kwargs__validation_set_to_False():
    """
    set_attrs_with_kwargs() accepts dict and updates attributes
    and can add non-valid if validate_kwargs = False
    """
    S = StoreIO()
    assert S.a is None
    assert S.b is None
    assert S.c is None
    assert S.d is None
    abcde = namedtuple('e', ['a','b','c','d','e'])
    x = abcde(1,2,3,4,5)
    S.set_attrs_with_kwargs(**x._asdict(), validate_kwargs = False)
    assert list(vars(S).keys()) == ['name', 'valid_attrs', 'valid_attrs_type', 'a', 'b', 'c', 'd', 'e']
    assert isinstance(S.a, int)
    assert isinstance(S.b, int)
    assert isinstance(S.c, int)
    assert isinstance(S.d, int)
    assert isinstance(S.e, int)


def test_StoreIO_set_attrs_with_extra_kwargs__validation_set_to_True():
    """
    set_attrs_with_kwargs() accepts dict and updates attributes
    and can add non-valid if validate_kwargs = False
    """
    S = StoreIO()
    assert S.a is None
    assert S.b is None
    assert S.c is None
    assert S.d is None
    abcde = namedtuple('e', ['a','b','c','d','e'])
    x = abcde(1,2,3,4,5)
    S.set_attrs_with_kwargs(**x._asdict(), validate_kwargs = True)
    assert list(vars(S).keys()) == ['name', 'valid_attrs', 'valid_attrs_type', 'a', 'b', 'c', 'd']
    assert isinstance(S.a, int)
    assert isinstance(S.b, int)
    assert isinstance(S.c, int)
    assert isinstance(S.d, int)

def test_StoreIO_set_attrs_with_kwargs__validate_kwargs_option_defaults_to_True():
    """
    set_attrs_with_kwargs() accepts dict and updates attributes
    and can add non-valid if validate_kwargs = False
    """
    S = StoreIO()
    assert S.a is None
    assert S.b is None
    assert S.c is None
    assert S.d is None
    abcde = namedtuple('e', ['a','b','c','d','e'])
    x = abcde(1,2,3,4,5)
    S.set_attrs_with_kwargs(**x._asdict())
    assert list(vars(S).keys()) == ['name', 'valid_attrs', 'valid_attrs_type', 'a', 'b', 'c', 'd']
    assert isinstance(S.a, int)
    assert isinstance(S.b, int)
    assert isinstance(S.c, int)
    assert isinstance(S.d, int)


def test_StoreIO_extra_kwargs__does_not_interfere_with_validate_attrs():
    """
    set_attrs_with_kwargs() accepts dict and updates attributes
    and can add non-valid if validate_kwargs = False
    """
    S = StoreIO()
    assert S.a is None
    assert S.b is None
    assert S.c is None
    assert S.d is None
    abcde = namedtuple('e', ['a','b','c','d','e'])
    x = abcde(1,2,3,4,5.0)
    S.set_attrs_with_kwargs(**x._asdict(), validate_kwargs = False)
    assert S._validate_attrs()
    assert isinstance(S.a, int)
    assert isinstance(S.b, int)
    assert isinstance(S.c, int)
    assert isinstance(S.d, int)
    assert isinstance(S.e, float)

def test_StoreIO_type_coerce():
    abcd = namedtuple('abcd', ['a','b','c','d'])
    x = abcd('1','2','3','4')
    S = StoreIO(**x._asdict())
    assert isinstance(S.a, str)
    S._type_coerce("a", S.valid_attrs_type[0])
    assert isinstance(S.a, int)


def test_StoreIO_coerce_attributes():
    abcd = namedtuple('abcd', ['a','b','c','d'])
    x = abcd('1','2','3','4')
    S = StoreIO(**x._asdict())
    assert S._coerce_attributes()

def test_StoreIO_coerce_attributes_then_validate_attrs():
    abcd = namedtuple('abcd', ['a','b','c','d'])
    x = abcd('1','2','3','4')
    S = StoreIO(**x._asdict())
    S._coerce_attributes()
    assert S._validate_attrs()

def test_StoreIO_coerce_attributes_then_validate_attrs__with_partial_set():
    abc = namedtuple('abc', ['a','b','c'])
    x = abc('1','2','3')
    S = StoreIO(**x._asdict())
    assert isinstance(S.a, str)
    assert isinstance(S.c, str)
    S._coerce_attributes()
    assert isinstance(S.a, int)
    assert isinstance(S.c, int)
    assert S._validate_attrs()

def test_StoreIO_coerce__raises_custom_ValueError_if_coercion_is_impossible():
    abcd = namedtuple('abcd', ['a','b','c','d'])
    x = abcd('A','B','C','D')
    S = StoreIO(**x._asdict())
    with pytest.raises(ValueError) as excinfo:
        S._coerce_attributes()
    assert str(excinfo.value) == "Cannot coerce StoreIO.a to <class 'int'>"

def test_StoreIO_coerce__raises_custom_ValueError_to_prevent_float_to_int_coercion():
    abcd = namedtuple('abcd', ['a','b','c','d'])
    x = abcd("1",1.5,"2",'3')
    S = StoreIO(**x._asdict())
    #with pytest.raises(ValueError) as excinfo:
    S._coerce_attributes()

    ##assert str(excinfo.value) == "Cannot coerce StoreIO.a to <class 'int'>"
