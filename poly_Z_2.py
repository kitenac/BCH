
'''
        List of crusual functions with poly_Z_2:

            0) there`s  +, -, *, /, % for polynoms
                also next() - itterator for polynoms

            1) Z2_prime_factors() - find all prime factors of poly (but dont count their deg: ex: x^2 devides F(x) => x will be appended, due x^2 - isn`t prime) <=> canonical factor without degrees
            2) get_primitive(t)   - gives primitive poly of deg = t that can permutate GF(2^t)/{0} and returns GF(2^t)/{0}`s int koeffs for it`s polynoms  
            3) on_val(g_x)        - returns value of f on g(x) <=> f(g(x))

'''




import copy              # for "phisical"(not by link) copping of objects   



# polynoms as set of bits
class poly_Z_2:


    # Problem:
    #   U can`t make multiple same-named functions in python. So the solution for polymorphism - named params or unnecessory(optional) params
    # Solution:
    #   - using "named" arguements except "positional":
    #   all parameters after * can be filled ONLY by their names. ex  __init__(koeffs = Idk, F_x = g.f_x)  - see https://tproger.ru/translations/asterisks-in-python-what-they-are-and-how-to-use-them/  
    # Feature: 
    #   "koeffs = None" <=> we can skip parametr, when calling function(due it have default init)  
    
    def __init__(self, *, koeffs = None, F = None, f_x:int = -1, str_koeffs = "" ):            # init poly by list of coefficents or F_x or ....
        
        if F != None:                           # by ecxisting polynom
            self.f_x = F.f_x
            self.deg = F.deg


        elif koeffs != None:                    # by list of koeffs
            self.f_x = 0
            self.deg = 1

            for i in range(0, len(koeffs)):
                self.f_x += koeffs[i] * 2**i

            self.update_deg()   # count degree

        elif f_x > -1:                          # by integer, which in binary represents polynom
            self.f_x = f_x
            self.deg = 0        # just init value
            self.update_deg()   # count degree

        elif str_koeffs != "":                  # by String of zeroes and ones
            lst = [int(k) for k in str_koeffs]
            self.__init__(koeffs=lst)

        else:                                   # default constructor
            self.deg = 0
            self.f_x = 0



    def __add__(self, G):
        S = poly_Z_2()
        
        S.f_x = self.f_x ^ G.f_x 
        S.update_deg()
        return S                #  "+" in GF_2 <=> XOR


    def update_deg(self):
        t_x = self.f_x 
        t_deg = 0

        while t_x > 1:      # when t_x = 1 or 0  - degree is 0
            t_x   //= 2 
            t_deg += 1
        
        self.deg = t_deg

        
        
    def __mul__(self, G):
        
        mul = poly_Z_2()
        for deg in range(0, self.deg + 1):
            
            if (self.f_x ^ 2**deg) < self.f_x:          # check if x^i excists (koeff by x^i = 1 )
                mul.f_x ^= G.f_x << deg                 # XOR is same as "+" in GF2

        mul.update_deg()    # due GF2 some powers could kill themself( 2 = 0 in GF2)
        return mul





    def __mod__(self, G):
        
        # Подвинули битово влево на разницу степеней и поксорили
        # step of mod-deviding
        def sub_mod(Res : poly_Z_2, f_mod):
            diff = Res.deg - f_mod.deg         # find optimal degree of x

            q_i  = f_mod.f_x << diff           # get f_mod * x^diff <=> rotate left "diff"-count bits  
            Res.f_x ^= q_i                     # XOR is "-" in GF_2          

            Res.update_deg()


        Res = poly_Z_2(F=self)                 # modulo-Remainder 

        if self.deg < G.deg:
            return Res
        else:
            while Res.deg >= G.deg:
                sub_mod(Res, G) 
        

        return Res


    def __truediv__(self, G):
        
        # Подвинули битово влево на разницу степеней и поксорили
        # step of mod-deviding
        def sub_div(Res : poly_Z_2, f_mod):
            diff = Res.deg - f_mod.deg         # find optimal degree of x
            q.f_x ^= 2**diff 

            q_i  = f_mod.f_x << diff           # get f_mod * x^diff <=> rotate left "diff"-count bits  
            Res.f_x ^= q_i                     # XOR is "-" in GF_2          
            Res.update_deg()

        Res = poly_Z_2(F=self)                 # modulo-Remainder 
        q = poly_Z_2()
        
        while Res.deg >= G.deg:
            sub_div(Res, G) 
            
        q.update_deg()
        return q







    # Phisical coppy(not by link) of poly. 
    # !!! Use it except "=" 
    def coppy(self):
        return poly_Z_2(F=self)
    
    def get_deg(self):
        return self.deg
    
    def get_f_x(self):
        return self.f_x


    def __str__(self):
        
        if self.deg == 0:
            return f"{self.f_x % 2}"

        f_x_str = ""                                        
        koefs = bin(self.f_x).lstrip("0b")                  # string with binary representation of f_x(stores as number) without format-string: "0b"
        k_deg = self.deg
    

        for koef in koefs:
            if koef != "0":  f_x_str += f"x^{k_deg} + "
            k_deg -= 1

        # pretty output for 1-th and 0-th powers
        f_x_str = f_x_str.replace("x^1 ", "x ")
        f_x_str = f_x_str.replace("x^0 ", "1 ")

        return f_x_str.rstrip(" + ") 


    def to_list(self):
        
        if self.deg == 0:
            return [self.f_x % 2]

        koefs = bin(self.f_x).lstrip("0b")                  # string with binary representation of f_x(stores as number) without format-string: "0b"
        f_x_lst = [int(k) for k in koefs]

        return f_x_lst[ : :-1]                              # return inverse ordered list               
        

    # proto iterator for polynoms
    def next(self):
        F_i = copy.deepcopy(self)
        F_i.f_x += 1                          # itterate over polynoms:     x, x+1, x^2, x^2 + 1, x^2 + x + 1, x^3, .... 
        if F_i.f_x % 2**(F_i.deg + 1) == 0:   # check if degree has grown (phisicaly)
            F_i.deg += 1                      # renew degree

        return F_i

    # tells if polynom is prime
    def is_prime(self):
        F_i = poly_Z_2(koeffs=[0, 1])             # F_i = x

        while F_i.deg <= self.deg//2:
            if  (self % F_i).f_x == 0:
                return False

            F_i = F_i.next()
        return True


    def uniqals(lst):
        uniq = [lst[0]]
        for el in lst:
            if el.f_x not in [u.f_x for u in uniq]:   uniq.append(el)
        return uniq


    
    def sub_factor(F, F_i):

        prime_factors = []
        
        while F_i.deg <= F.deg - 1:

            if F_i.is_prime() and (F % F_i).f_x == 0:
                
                prime_factors.append( copy.deepcopy(F_i) )                 # adding phisical coppy of F_i, bc object F_i will be modified later <=> F_i in list will change too 
                sub = poly_Z_2.sub_factor( (F / F_i), F_i.next())          # now factorizing subfactor
                
                prime_factors.extend(sub)
                prime_factors = poly_Z_2.uniqals(prime_factors)            # returns only uniqual factors  
                return prime_factors

            F_i = F_i.next()


        if len(prime_factors) == 0:
            prime_factors.append(F)
            return  prime_factors


    # returns only prime deviders of poly. ex: Z2_factorize(x^4 + x^2) = x, x^2 + 1 - not x^2, x^2 + 1, bc x^2 - not prime
    def Z2_prime_factors(self):
        F_i = poly_Z_2(koeffs=[0, 1])             # F_i = x - value to start reqursion

        return poly_Z_2.sub_factor(self, F_i)



    def pow(self, k: int):
        if k == 0:   return poly_Z_2(koeffs=[1])
        tmp = copy.deepcopy(self)
        for _ in range(1, k):
            tmp *= self
        return tmp

    def on_value(self, val):
        koefs = self.to_list()
        res = poly_Z_2()
        for i in range (0, len(koefs)):
            if koefs[i] != 0:   res += val.pow(i)

        return res   

    # if deg = 13, 12 - 2 sec / 0.7 sec, if deg <= 11 - just a moment. Otherwise - too long
    def get_primitive(deg: int):
        #print(f"Let`s find primitive poly of degree = {deg}" )

        f_i = poly_Z_2(f_x=2**deg)
        alpha = poly_Z_2(f_x=2)       # x is always primitive
        while True:
            f_i = f_i.next() 
            if f_i.is_prime():
                #print(f"Check GF/{f_i}, f_i = {f_i}:")
                test_GF = [1]                       # may be GF
                el = poly_Z_2(f_x=1)                # x^0 = 1 
                while len(test_GF) < 2**deg - 1:    # while GF not fully filled
                    el = (el * alpha) % f_i 
                    #print(f"({el} * {alpha}) % f_i  = {(el * alpha) % f_i}")
                    if el.f_x not in test_GF:
                        test_GF.append(el.f_x)
                    else: 
                        #print(f" !!! element {(el * alpha) % f_i} is already in GF => go to next prime poly:\n")
                        break                       # found same element => it`s not GF => take next poly
                else:
                    #print(f"\n-Finaly we build GF/{f_i}")
                    return f_i, test_GF             # if test_GF is GF(consists of uniqual elements). ps test_GF 

            if f_i.deg > deg:
                print(f"Something went wrong(BUG) and there is no primitive poly(but must be) of deg = {deg}")    



# ---------------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    coefs = [1, 0, 1, 1, 1]
    koefs = [0, 0, 0, 1, 0, 0, 1, 1]

    F = poly_Z_2(koeffs = coefs)
    G = poly_Z_2(koeffs = koefs)

    print(f"F = {F},    G = {G} \n F+G = {(F + G)},\n G mod F = {(G % F)}, degree = {(G % F).deg} " )

    print("\nChecking f_x --> list:")
    print (f" {F}    as list: {F.to_list()} \n {G}   as list {G.to_list()}") 
    
    print(f"F*G = {F * G}")


    A = poly_Z_2(koeffs=[1, 1, 0, 1])
    B = poly_Z_2(koeffs=[1, 1, 0, 0])
    print(f"A = {A}, B = {B},\n A/B = {A / B},  A mod B = {A%B} ")

  
    print(f"\n\n Factorization of {G.next()}: ")
    for el in G.next().Z2_prime_factors():
        print(f"{el}", end = ", ")
    print("\n")



    #print(f"\n\nF(x) = {F},\nF(x+1) = {F.on_value( poly_Z_2(koeffs=[1, 1]) )} ")
    

    #print("Let`s find primitive polynom of deg = t.\n  ps this polynom will permutate GF(2^t)/\{0\}\n")

    #for t in range(2, 13):
    #    print(f"t = {t}: {poly_Z_2.get_primitive(t)[0]}")


    print(f"\ng(x) = {poly_Z_2.get_primitive(8)[0]}" )