"""Directed MPFR enclosures for mathematical normal CDF, never Excel answers.

Uses the installed MPFR shared library, not its source. No code is copied from
MPFR. This is a mathematical diagnostic, not an x87 emulator or a search rank.
Only the official MPFR interface/rounding contract and DLMF7.20.1 are used.
"""
import ctypes as C
import ctypes.util
from fractions import Fraction as R
import json
from pathlib import Path
import struct

from exact import BINARY64 as B, binary64, two
from check import require, digest, records, write_result


class Raw(C.Structure):
    _fields_ = [('precision', C.c_long), ('sign', C.c_int),
                ('exponent', C.c_long), ('limbs', C.POINTER(C.c_ulong))]


LIB = C.CDLL(ctypes.util.find_library('mpfr'))
PTR = C.POINTER(Raw)
for name, args, result in [
        ('init2', [PTR,C.c_long], None), ('clear',[PTR],None),
        ('set_d',[PTR,C.c_double,C.c_int],C.c_int),
        ('sqrt',[PTR,PTR,C.c_int],C.c_int), ('erfc',[PTR,PTR,C.c_int],C.c_int),
        ('log',[PTR,PTR,C.c_int],C.c_int),
        ('jn',[PTR,C.c_long,PTR,C.c_int],C.c_int),
        ('yn',[PTR,C.c_long,PTR,C.c_int],C.c_int),
        ('div',[PTR,PTR,PTR,C.c_int],C.c_int), ('mul',[PTR,PTR,PTR,C.c_int],C.c_int),
        ('sub',[PTR,PTR,PTR,C.c_int],C.c_int),
        ('get_d',[PTR,C.c_int],C.c_double),
        ('get_str',[C.c_void_p,C.POINTER(C.c_long),C.c_int,C.c_size_t,PTR,C.c_int],C.c_void_p),
        ('free_str',[C.c_void_p],None), ('get_version',[],C.c_char_p)]:
    f = getattr(LIB, 'mpfr_'+name)
    f.argtypes, f.restype = args, result


def rational(value):
    exponent = C.c_long()
    ptr = LIB.mpfr_get_str(None,C.byref(exponent),2,0,C.byref(value),0)
    try:
        digits = C.string_at(ptr).decode()
        require(digits.strip('-01') == '', 'nonfinite MPFR value')
        negative = digits.startswith('-')
        digits = digits.lstrip('-')
        return (-1 if negative else 1)*int(digits,2)*two(exponent.value-len(digits))
    finally:
        LIB.mpfr_free_str(ptr)


def log_bounds(bits, precision=256):
    """Directed MPFR bounds on log of an exact positive binary64 input."""
    value=struct.unpack('>d',int(bits).to_bytes(8,'big'))[0]
    require(value>0 and value<float('inf'),'positive finite logarithm input required')
    raw=[Raw(),Raw(),Raw()]
    for item in raw:LIB.mpfr_init2(C.byref(item),precision)
    try:
        LIB.mpfr_set_d(C.byref(raw[0]),value,0)
        LIB.mpfr_log(C.byref(raw[1]),C.byref(raw[0]),3)
        LIB.mpfr_log(C.byref(raw[2]),C.byref(raw[0]),2)
        return rational(raw[1]),rational(raw[2])
    finally:
        for item in raw:LIB.mpfr_clear(C.byref(item))


def bessel_bounds(function, bits, order, precision=256):
    """Directed real-function bounds through the official MPFR jn/yn API.

    Calls the installed shared library only; no MPFR implementation is read.
    These are mathematical references, never Excel answers or native x87 leaves.
    """
    require(function in ('BESSELJ','BESSELY') and 0 <= order <= 64, 'Bessel reference scope')
    value=struct.unpack('>d',int(bits,16).to_bytes(8,'big'))[0]
    values=[Raw(),Raw(),Raw()]
    for raw in values:LIB.mpfr_init2(C.byref(raw),precision)
    try:
        LIB.mpfr_set_d(C.byref(values[0]),value,0)
        operation=LIB.mpfr_jn if function=='BESSELJ' else LIB.mpfr_yn
        operation(C.byref(values[1]),order,C.byref(values[0]),3)
        operation(C.byref(values[2]),order,C.byref(values[0]),2)
        lo,hi=rational(values[1]),rational(values[2])
        require(lo<=hi,'reversed Bessel reference interval')
        return lo,hi
    finally:
        for raw in values:LIB.mpfr_clear(C.byref(raw))


def normal_bounds(bits, precision=256):
    """Propagate directed bounds through sqrt/div/erfc/half/complement."""
    bits = int(bits,16) if isinstance(bits,str) else bits
    raw = []
    def new(number=None):
        r = Raw();LIB.mpfr_init2(C.byref(r),precision);raw.append(r)
        if number is not None:
            LIB.mpfr_set_d(C.byref(r),number,0)
        return r
    def call(name, *args):
        r = new();getattr(LIB,'mpfr_'+name)(C.byref(r),*(C.byref(a) for a in args[:-1]),args[-1]);return r
    x = binary64(bits)
    require(abs(x)<=40, 'normal reference restricted to |z|<=40')
    try:
        argument = new(struct.unpack('>d',(bits&((1<<63)-1)).to_bytes(8,'big'))[0])
        two_value, half, one = new(2.0), new(0.5), new(1.0)
        sqrt_low, sqrt_high = call('sqrt',two_value,3), call('sqrt',two_value,2)
        t_low = call('div',argument,sqrt_high,3)
        t_high = call('div',argument,sqrt_low,2)
        lower = call('mul',call('erfc',t_high,3),half,3)
        upper = call('mul',call('erfc',t_low,2),half,2)
        if x>0:
            lower,upper = call('sub',one,upper,3),call('sub',one,lower,2)
        lo,hi = rational(lower),rational(upper)
        require(lo<=hi,'reversed directed interval')
        low_bits = struct.pack('>d',LIB.mpfr_get_d(C.byref(lower),0)).hex()
        high_bits = struct.pack('>d',LIB.mpfr_get_d(C.byref(upper),0)).hex()
        require(low_bits==high_bits,'insufficient mathematical reference precision')
        require(B.round(lo)==B.round(hi)==binary64(low_bits),'MPFR conversion disagrees with rational lattice')
        return dict(lower=[str(lo.numerator),str(lo.denominator)],upper=[str(hi.numerator),str(hi.denominator)],
                    rounded_binary64='0x'+low_bits,precision=precision)
    finally:
        for r in raw:LIB.mpfr_clear(C.byref(r))


def diagnostic(trace, destination):
    from collections import Counter
    destination = Path(destination)
    destination.mkdir(exist_ok=True)
    counts = Counter()
    with (destination/'math-cdf.jsonl').open('x') as out:
        for index,row in enumerate(records(trace)):
            a = row['args']
            flag = a[1] if isinstance(a[1],bool) else binary64(a[1])
            if not flag:
                counts['density_control']+=1;continue
            bits = a[0]
            if abs(binary64(bits))>40:
                counts['outside_reference_domain']+=1;continue
            ref = normal_bounds(bits)
            build = row['reading']['excel_build']
            match = ref['rounded_binary64'] == row['expected']
            counts[f'{build}:mathematical_rounding_equal' if match else f'{build}:mathematical_rounding_differs']+=1
            # Predeclared transport hypotheses; counts never select an algorithm.
            alternatives = {}
            value = struct.unpack('>d',bytes.fromhex(bits[2:]))[0]
            for digits in (15,16,17):
                literal = format(value,f'.{digits}g')
                parsed = '0x'+struct.pack('>d',float(literal)).hex()
                alternative = normal_bounds(parsed) if parsed.lower()!=bits.lower() else ref
                alternatives[str(digits)] = dict(literal=literal,parsed_by_python=parsed,
                    expected_equal_math=alternative['rounded_binary64']==row['expected'])
            out.write(json.dumps(dict(trace_index=index,id=row['id'],build=build,input_bits=bits,
                expected=row['expected'],reference=ref,diagnostic_decimal_transports=alternatives))+'\n')
    result = dict(counts=counts,mpfr_version=LIB.mpfr_get_version().decode(),
        reference_source_sha256=digest(Path(__file__)),trace_sha256=digest(Path(trace)),
        claim='Directed mathematical CDF bounds; deviations diagnose assumptions, not program admissibility or literal Excel inputs. Decimal alternatives use Python parsing, not measured Excel parsing.',
        new_excel_answers=0,algorithm_selected=False)
    write_result(destination/'summary.json',result)


if __name__=='__main__':
    import sys
    diagnostic(*sys.argv[1:])
