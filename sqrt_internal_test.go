package bigfloat

import (
	"math/big"
	"math/rand"
	"testing"
)

// TestSqrtNewtonAccuracy holds the root to 3 ulps of w, on which both Sqrt's
// halfway test and the AGM loops rely.
func TestSqrtNewtonAccuracy(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	for _, w := range []uint{1, 53, 65, 88, 89, 101, 117, 164, 228, 1000, 6664, 13364, 26692, 54000} {
		for range 40 {
			xp := []uint{w, 2*w + 10, 53, 4 * w}[rng.Intn(4)]
			m := new(big.Int).Rand(rng, new(big.Int).Lsh(big.NewInt(1), xp))
			x := new(big.Float).SetPrec(xp).SetInt(m.Add(m, big.NewInt(1)))
			x.SetMantExp(x, rng.Intn(2000)-1000-int(xp))
			got := sqrtNewton(x, w)
			want := sqrtExact(new(big.Float).SetPrec(w+40), x)
			d := new(big.Float).SetPrec(w+60).Sub(got, want)
			ulp := new(big.Float).SetMantExp(big.NewFloat(1), want.MantExp(nil)-int(w))
			if r, _ := d.Abs(d).Quo(d, ulp).Float64(); r > 3 {
				t.Errorf("w = %d: sqrtNewton(%.20g) is %.3g ulps off", w, x, r)
			}
		}
	}
}
