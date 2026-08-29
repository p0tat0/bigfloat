package bigfloat

import "math/big"

// Pow returns a big.Float representation of zʷ. Precision is the same
// as the one of the first argument. The function panics when z is
// negative and w is not an integer.
func Pow(z *big.Float, w *big.Float) *big.Float {

	if z.Sign() < 0 {
		// For negative base, handle integer exponents: (-z)^n = (-1)^n * z^n
		if w.IsInt() {
			// Compute |z|^w
			abs := new(big.Float).Abs(z)
			result := Pow(abs, w)
			// If exponent is odd, negate the result
			wInt, _ := w.Int(nil)
			if wInt.Bit(0) == 1 { // odd
				result.Neg(result)
			}
			return result
		}
		panic("Pow: negative base with non-integer exponent")
	}

	// Pow(z, 0) = 1.0
	if w.Sign() == 0 {
		return big.NewFloat(1).SetPrec(z.Prec())
	}

	// Pow(z, 1) = z
	// Pow(+Inf, n) = +Inf
	if w.Cmp(big.NewFloat(1)) == 0 || z.IsInf() {
		return new(big.Float).Copy(z)
	}

	// Pow(z, -w) = 1 / Pow(z, w)
	if w.Sign() < 0 {
		x := new(big.Float)
		zExt := new(big.Float).Copy(z).SetPrec(z.Prec() + 64)
		wNeg := new(big.Float).Neg(w)
		return x.Quo(big.NewFloat(1), Pow(zExt, wNeg)).SetPrec(z.Prec())
	}

	// compute w**z as exp(z log(w))
	x := new(big.Float).SetPrec(z.Prec() + 64)
	logZ := Log(new(big.Float).Copy(z).SetPrec(z.Prec() + 64))
	x.Mul(w, logZ)
	x = Exp(x)
	return x.SetPrec(z.Prec())

}
