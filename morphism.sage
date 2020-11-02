
_.<X> = GF(2)[]
K.<x> = GF(2^8, modulus= X^8 +X^4+X^3+X+1)


F = GF(2^32, 'y', modulus = X^32 + X^22 + X^2 + X^1 + 1)
print(F, F.modulus())
print(F.gens())
H = Hom(K, F)
morph = H.list()[0]
g = morph.im_gens()[0]
print(morph)
#a = K.random_element()
#print(g^7 + g^5+g^3+1 == morph(a))
