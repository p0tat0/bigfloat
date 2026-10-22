package bigfloat_test

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/p0tat0/bigfloat"
)

var roundingModes = []big.RoundingMode{
	big.ToNearestEven, big.ToNearestAway, big.ToZero,
	big.AwayFromZero, big.ToNegativeInf, big.ToPositiveInf,
}

// checkRounded fails unless r is √x correctly rounded to prec bits in mode,
// certified by exact squares of r's rounding interval.
func checkRounded(t *testing.T, x, r *big.Float, prec uint, mode big.RoundingMode) {
	t.Helper()
	mant := new(big.Float)
	e := r.MantExp(mant)
	ulp := new(big.Float).SetMantExp(big.NewFloat(1), e-int(prec))
	below := new(big.Float).Set(ulp)
	if mant.Cmp(big.NewFloat(0.5)) == 0 {
		below.SetMantExp(below, -1)
	}
	sq := func(v *big.Float) *big.Float {
		return new(big.Float).SetPrec(2*prec+8).Mul(v, v)
	}
	at := func(k int64, step *big.Float) *big.Float {
		d := new(big.Float).SetPrec(prec+8).Mul(step, big.NewFloat(float64(k)))
		d.SetMantExp(d, -1)
		return new(big.Float).SetPrec(prec+8).Add(r, d)
	}
	var lo, hi *big.Float
	var loOpen, hiOpen bool
	switch mode {
	case big.ToNearestEven, big.ToNearestAway:
		lo, hi = sq(at(-1, below)), sq(at(1, ulp))
	case big.ToZero, big.ToNegativeInf:
		lo, hi, hiOpen = sq(r), sq(at(2, ulp)), true
	default:
		lo, hi, loOpen = sq(at(-2, below)), sq(r), true
	}
	if c := lo.Cmp(x); c > 0 || (loOpen && c == 0) {
		t.Errorf("prec %d mode %v: Sqrt(%.20g) = %.20g is too large", prec, mode, x, r)
	}
	if c := x.Cmp(hi); c > 0 || (hiOpen && c == 0) {
		t.Errorf("prec %d mode %v: Sqrt(%.20g) = %.20g is too small", prec, mode, x, r)
	}
}

// TestSqrtCorrectlyRounded covers precisions where math/big's root is a few
// to 2^127 ulps off: prec + 32 just below 53·2^k.
func TestSqrtCorrectlyRounded(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for _, prec := range []uint{24, 53, 64, 200, 1000, 13300, 26628} {
		inputs := []*big.Float{}
		for _, i := range []int{70, 102, 108, 189} {
			x := new(big.Float).SetPrec(prec).SetInt64(int64(i))
			x.Quo(x, new(big.Float).SetPrec(prec).SetInt64(201))
			inputs = append(inputs, x.SetMantExp(x, -(i%7)))
		}
		for range 12 {
			xp := []uint{prec, 2*prec + 10, 53}[rng.Intn(3)]
			m := new(big.Int).Rand(rng, new(big.Int).Lsh(big.NewInt(1), xp))
			x := new(big.Float).SetPrec(xp).SetInt(m.Add(m, big.NewInt(1)))
			inputs = append(inputs, x.SetMantExp(x, rng.Intn(400)-200-int(xp)))
		}
		for _, x := range inputs {
			for _, mode := range roundingModes {
				z := new(big.Float).SetPrec(prec).SetMode(mode)
				checkRounded(t, x, bigfloat.Sqrt(z, x), prec, mode)
			}
		}
	}
}

// TestSqrtHalfway: for an odd (prec+1)-bit r, √(r²) is halfway between two
// prec-bit values, and √(r² ± 1) is within 2^-(2 prec) of halfway, closer
// than the refined root resolves.
func TestSqrtHalfway(t *testing.T) {
	const prec = 100
	rng := rand.New(rand.NewSource(2))
	for range 20 {
		r := new(big.Int).Rand(rng, new(big.Int).Lsh(big.NewInt(1), prec))
		r.SetBit(r, prec, 1).SetBit(r, 0, 1)
		down := new(big.Float).SetPrec(prec).SetMode(big.ToZero).SetInt(r)
		up := new(big.Float).SetPrec(prec).SetMode(big.AwayFromZero).SetInt(r)
		for _, mode := range []big.RoundingMode{big.ToNearestEven, big.ToNearestAway} {
			for d, want := range map[int64]*big.Float{
				-1: down,
				0:  new(big.Float).SetPrec(prec).SetMode(mode).SetInt(r),
				1:  up,
			} {
				sq := new(big.Int).Mul(r, r)
				x := new(big.Float).SetPrec(2*prec + 2).SetInt(sq.Add(sq, big.NewInt(d)))
				got := bigfloat.Sqrt(new(big.Float).SetPrec(prec).SetMode(mode), x)
				if got.Cmp(want) != 0 {
					t.Errorf("mode %v: Sqrt(%v² %+d) = %v, want %v", mode, r, d, got, want)
				}
			}
		}
	}
}

func TestSqrtSpecialValues(t *testing.T) {
	for _, x := range []*big.Float{
		new(big.Float), new(big.Float).Neg(new(big.Float)), new(big.Float).SetInf(false),
		big.NewFloat(4), big.NewFloat(0.25),
	} {
		got := bigfloat.Sqrt(new(big.Float), x)
		want := new(big.Float).Sqrt(x)
		if got.Cmp(want) != 0 || got.Signbit() != want.Signbit() || got.Prec() != x.Prec() {
			t.Errorf("Sqrt(%v) = %v (prec %d), want %v (prec %d)", x, got, got.Prec(), want, x.Prec())
		}
	}
	x := big.NewFloat(2)
	if got := bigfloat.Sqrt(x, x); got.Cmp(big.NewFloat(1.4142135623730951)) != 0 {
		t.Errorf("Sqrt(x, x) with x = 2: %v", got)
	}
	defer func() {
		if recover() == nil {
			t.Error("Sqrt(-1) did not panic")
		}
	}()
	bigfloat.Sqrt(new(big.Float), big.NewFloat(-1))
}
