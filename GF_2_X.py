
import copy 
from poly_Z_2 import poly_Z_2

import random


'''
  Reduced pole: GF_2_X\{0} - there is no  0 element( useless in BCH and can`t be permute from alpha by def)
   !!! Also GF_2_X\{0} - is multiplicative group of pole

    Some crusual functions:
        1) get_ord(), get_alpha_i() - get order of some poly in some GF
        2) mul(), mod(), dev()      - *, %, / - but cheaper - in GF
        3) get_g_x()                - get generator poly: g(alpha) = 0 for each apha in GF
'''


class GF2_x:

    def __init__(self, P: poly_Z_2, *, Table: list = None):
        self.alpha = poly_Z_2(koeffs = [0, 1])       # alpha = x is always primitive, so let it be built-in  
        self.p_x = P.coppy()

        # Table of int-represented polynoms - useful with poly_Z_2.get_primitive.  This case - much cheapper
        if Table != None:
            self.table = [poly_Z_2(f_x=el) for el in Table]

        # counting each el of table
        else:
            self.table = [ poly_Z_2(koeffs = [1]) ]      # 0-th power of alpha = 1 is always in table(GF_X) 
            for i in range(1, 2**P.get_deg() - 1):      
                self.table.append(  (self.table[i-1] * self.alpha) % P)  # counting each alpha^i = alpha^i-1 * alpha mod p_x
        
        




    def get_alpha_i(self, i: int):
        return copy.deepcopy( self.table[i % len(self.table)] )

    def get_ord(self, v : poly_Z_2):
        i = 0
        for f in self.table:
            if f.f_x == v.f_x:  return i
            i += 1    
        print(f"\n ==================== !!Alert!! ===================\n  {self}\n !!Alert!! There is no ({v}) in GF(watch table above)")    

    
    def get_table(self):
        return self.table 


        # F^-1: F * F^-1 = 1
    def get_opposite(self, F:poly_Z_2):
        ord = self.get_ord(F)
        return copy.deepcopy( self.get_alpha_i(len(self.table) - ord) )


    # Cheap multiplying in GF_2^t using table:
    def mul(self, f: poly_Z_2, g: poly_Z_2):
        if f.f_x == 0 or g.f_x == 0: 
            return poly_Z_2(f_x=0)

        i, j = self.get_ord(f),  self.get_ord(g) 
        return copy.deepcopy( self.get_alpha_i( (i+j) % len(self.table)) )


    def pow(self, F: poly_Z_2, k: int):
        if k == 0:   
            return poly_Z_2(koeffs=[1])

        ord = self.get_ord(F)
        return self.get_alpha_i( (ord * k) % len(self.table)) 


    def devide(self, F: poly_Z_2, G: poly_Z_2):
        if F.deg < G.deg:
            return poly_Z_2()
        return copy.deepcopy( self.get_alpha_i( (F.deg - G.deg) % len(self.table) ))



    def mod(self, F: poly_Z_2, G: poly_Z_2):
        return F - self.mul(G, self.devide(F, G))    # F - (F//G)*G . F//G - is   


    def on_val(self, F, val):
        koefs = F.to_list()
        res = poly_Z_2()
        for i in range (0, len(koefs)):
            if koefs[i] != 0:   res += self.pow(val, i)

        return res % self.p_x  


    # seek for minimal polyom g_x: g(alpha) == 0 <=> each alpha devides g(x)
    def get_g_x(self, t: int):
        devs = []               

        P = poly_Z_2(f_x = 2**(2**self.p_x.deg) + 2)        # f(x) = x^(2^m) + x - this poly always has deviders for g(x) in GF(2^m) 
        P_devs = P.Z2_prime_factors()

        delta = 2*t + 1                                     # minimal code distanse

        for g_root in self.table[1:delta]:
            for poly in P_devs:
                if self.on_val(poly, g_root).get_f_x() == 0:
                    devs.append(poly)
                    #print(f"   for α^{self.get_ord(g_root)} root = {poly}")
                    break

        #print("\nNow we reduce multiple deviders and get g(x):")
        devs = poly_Z_2.uniqals(devs)                       # Оставляем только уникальные делители <=> НОК

        g_x = poly_Z_2(f_x=1)
        for d in devs:      g_x *= d
        

        return g_x




    def __str__(self):
        deg = 0
        s = ""

        s +=  f"\n GF_2^{self.p_x.deg}, p(x) = {self.p_x}  |GF_2^{self.p_x.deg}\\0|  = {len(self.table)}\n\n"

        for alpha in self.table:
            s += f" α^{deg} = {alpha} \n"
            deg += 1
        
        return s
        

    








if __name__ == "__main__":
    GF2_3 = GF2_x(poly_Z_2(koeffs=[1, 1, 0, 1]))          # p(x) = x^3 + x + 1 - primitive poly => can permutate all GF`s element as x^i mod p(x) 
    GF2_4 = GF2_x(poly_Z_2(koeffs=[1, 1, 0, 0, 1]))       # x^4 + x + 1 - primitive poly
    GF2_5 = GF2_x(poly_Z_2(koeffs=[1, 0, 1, 0, 0, 1]))       # x^5 + x^2 + 1 - primitive poly


    print(f"\n{GF2_3} \n{GF2_4} \n{GF2_5}")

    print(f" alpha`s power to get x^2:  {GF2_3.get_ord(poly_Z_2(koeffs=[0,0,1]))}" )


    f = poly_Z_2(koeffs=[0, 1, 0])
    g = poly_Z_2(koeffs=[0, 1, 1])
    print(f" {f} * {g} = {f*g}     in GF_2, \n {f} * {g} = {GF2_3.mul(f, g)}    in GF_2/{GF2_3.p_x} \n {f} * {g} = {GF2_4.mul(f, g)}    in GF_2/{GF2_4.p_x}  ")


    
    print("\n\nLet`s find p(x) for BCH-coder: t = 2, GF = GF(2^4): ")
    print(GF2_4)
    G = GF2_4.get_g_x(2)    
    print(f"g(x) = {G}")

    #print(f"({g}) * ({GF2_3.get_opposite(g)}) = {g * GF2_3.get_opposite(g)}")






