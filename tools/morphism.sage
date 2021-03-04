
R.<X> = GF(2)[]
K.<x> = GF(2^8, modulus= X^8 +X^4+X^3+X+1)


F = GF(2^16, 'y', modulus = X^16 + X^5 + X^3 + X + 1)
#F = GF(2^32, 'y', modulus = X^32 + X^7 + X^3 + X^2 + 1)
#F = GF(2^40, 'y', modulus = X^40 + X^5 + X^4 + X^3 + 1)
#F = GF(2^48, 'y', modulus = X^48 + X^5 + X^3 + X^2 + 1)

print(F, F.modulus())
print(F.gens())
H = Hom(K, F)
morph = H.list()[0]
g = morph.im_gens()[0]
print(morph)
