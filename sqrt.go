package bigfloat

import (
	"math"
	"math/big"
)

// sqrtGuard is the precision Sqrt works at above z's. Its Newton iteration
// leaves the root within 3 ulps of √x there.
const sqrtGuard = 64

// Sqrt sets z to the square root of x, rounded to z's precision in z's
// rounding mode, and returns z. If z's precision is 0, it is changed to x's
// before the operation. Sqrt panics if x < 0.
//
// Unlike z.Sqrt(x), the result is correctly rounded. math/big's Newton steps
// double the precision from 53 bits and lose about a bit per step, so when
// prec + 32 lies just below 53·2^k up to 2% of the bits come out wrong.
func Sqrt(z, x *big.Float) *big.Float {
	if z.Prec() == 0 {
		z.SetPrec(x.Prec())
	}
	if x.Sign() < 0 {
		panic("Sqrt: argument is negative")
	}
	if x.Sign() == 0 || x.IsInf() {
		return z.Set(x)
	}
	prec := z.Prec()
	y := sqrtNewton(x, prec+sqrtGuard)
	if m := z.Mode(); (m != big.ToNearestEven && m != big.ToNearestAway) || nearHalfway(y, prec) {
		return sqrtExact(z, x)
	}
	return z.Set(y)
}

// three is read only.
var three = big.NewFloat(3)

// sqrtNewton is √x to w bits, within 3 ulps: 1/√m for the mantissa m in
// [1/2, 2) by t ← t(3 − m t²)/2 from a float64 seed, each step at a precision
// its predecessor's accuracy covers with 8 bits to spare, then m t. The
// seed's 52 bits carry the first step to 100.
func sqrtNewton(x *big.Float, w uint) *big.Float {
	m := new(big.Float)
	e := x.MantExp(m)
	if e%2 != 0 {
		m.SetMantExp(m, 1)
		e--
	}
	m.SetPrec(w + 8)
	var buf [64]uint
	steps := append(buf[:0], w)
	for p := w; p > 100; {
		p = p/2 + 8
		steps = append(steps, p)
	}
	f, _ := m.Float64()
	t := new(big.Float).SetFloat64(1 / math.Sqrt(f))
	u := new(big.Float)
	for i := len(steps) - 1; i >= 0; i-- {
		p := steps[i]
		t.SetPrec(p)
		u.SetPrec(p).Mul(t, t)
		u.Mul(u, m)
		u.Sub(three, u)
		t.Mul(t, u)
		t.SetMantExp(t, -1)
	}
	t.SetPrec(w).Mul(t, m)
	return t.SetMantExp(t, e/2)
}

// nearHalfway reports whether y, held to prec + sqrtGuard bits, lies within
// 16 of its ulps of halfway between two prec-bit neighbours, where its error
// could decide the rounding.
func nearHalfway(y *big.Float, prec uint) bool {
	c := new(big.Float).SetPrec(prec).SetMode(big.ToZero).Set(y)
	d := new(big.Float).SetPrec(sqrtGuard+1).Sub(y, c)
	e := c.MantExp(nil) - int(prec)
	half := new(big.Float).SetMantExp(big.NewFloat(1), e-1)
	d.Sub(d, half).Abs(d)
	return d.Sign() == 0 || d.MantExp(nil) <= e-int(sqrtGuard)+4
}

// sqrtExact sets z to √x from the integer square root of x's mantissa scaled
// to at least 2·prec + 4 bits, with a sticky bit for an inexact root, so the
// one rounding is z's.
func sqrtExact(z, x *big.Float) *big.Float {
	prec := z.Prec()
	mant := new(big.Float)
	e := x.MantExp(mant)
	s := int(max(x.Prec(), 2*prec+4))
	if (e-s)%2 != 0 {
		s++
	}
	m, _ := mant.SetMantExp(mant, s).Int(nil)
	r := new(big.Int).Sqrt(m)
	r.Lsh(r, 1)
	if new(big.Int).Mul(r, r).Cmp(m.Lsh(m, 2)) != 0 {
		r.SetBit(r, 0, 1)
	}
	z.SetInt(r)
	return z.SetMantExp(z, (e-s)/2-1)
}
