_.<X> = GF(2)[]
K.<x> = GF(2^8, modulus= X^8 +X^4+X^3+X+1)
F = GF(2^64, 'y')
H = Hom(K, F)
morph = H.list()[0]
g = morph.im_gens()[0]
a = K.random_element()
print(g^7 + g^5+g^3+1 == morph(a))
