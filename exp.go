package bigfloat

import (
	"math"
	"math/big"
)

// smallestNormalFloat64 is the smallest positive normal float64; below
// it math.Exp returns subnormals, whose precision drops from 53 bits to
// as few as one, and 0 below 2⁻¹⁰⁷⁴.
const smallestNormalFloat64 = 0x1p-1022

// Exp returns a big.Float representation of eᶻ, aka exp(z). Precision
// is the same as the one of the argument. Returns +Inf when z = +Inf,
// and 0 when z = -Inf.
func Exp(z *big.Float) *big.Float {

	// exp(0) == 1
	if z.Sign() == 0 {
		return big.NewFloat(1).SetPrec(z.Prec())
	}

	// Exp(+Inf) = +Inf
	if z.IsInf() && z.Sign() > 0 {
		return big.NewFloat(math.Inf(+1)).SetPrec(z.Prec())
	}

	// Exp(-Inf) = 0
	if z.IsInf() && z.Sign() < 0 {
		return big.NewFloat(0).SetPrec(z.Prec())
	}

	// |z| >= 2^31 is far past MaxExp·ln 2 and (MinExp-1)·ln 2; below
	// that, the halving recursion is shallow and decides the edge.
	if z.MantExp(nil) > 31 {
		if z.Sign() > 0 {
			return new(big.Float).SetPrec(z.Prec()).SetInf(false)
		}
		return new(big.Float).SetPrec(z.Prec())
	}

	guess := new(big.Float)

	// try to get initial estimate using IEEE-754 math
	zf, _ := z.Float64()
	if zfs := math.Exp(zf); zfs == math.Inf(+1) || zfs < smallestNormalFloat64 {
		// too big or too small for IEEE-754 math (overflow, underflow,
		// or a subnormal result carrying fewer than 53 correct bits —
		// newton below assumes a full-precision seed and would trust
		// it), perform argument reduction using
		//     e^{2z} = (e^z)²
		halfZ := new(big.Float).Mul(z, big.NewFloat(0.5))
		halfExp := Exp(halfZ.SetPrec(z.Prec() + 64))
		return new(big.Float).Mul(halfExp, halfExp).SetPrec(z.Prec())
	} else {
		// we got a nice IEEE-754 estimate
		guess.SetFloat64(zfs)
	}

	// f(t)/f'(t) = t*(log(t) - z)
	f := func(t *big.Float) *big.Float {
		x := new(big.Float)
		x.Sub(Log(t), z)
		return x.Mul(x, t)
	}

	x := newton(f, guess, z.Prec())

	return x
}
