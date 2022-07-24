

from copy import deepcopy

from typing import List
from poly_Z_2 import poly_Z_2
from GF_2_X import GF2_x
import random   



# Даты сдачи 27, 30, 31, (Г-401) 3,4 пары;  1, и не факт 2, 3


'''
 Polynoms(y) that has poly_Z_2 as koefficents by degrees of y. also I assume poly_Y`s koefs to bellong GF for easier multiplying(by primitive element powers) 

'''
class poly_Y:
        def __init__(self, GF: GF2_x, *, koeffs: List[poly_Z_2] = [], koeffs_alpha: List[int] = []):

            self.GF = GF
            
            # koeffs must be without extra koeffs: ex:  x*y^2 => koeffs = [0, 0, x],  not [0, 0, x, 0]
            if koeffs != []:
                self.koeffs = deepcopy(koeffs) 
                self.deg = len(koeffs) - 1
            
            elif koeffs_alpha != []:
                self.koeffs = [ GF.get_alpha_i(k) for k in koeffs_alpha ]       # fill koeffs of polynom(y) due GF table and given koeffs
                self.deg = len(koeffs_alpha) - 1

            # Creating emty poly_Y that can be filled appending elements
            else:
                self.koeffs = []    # nothing. even not a zero
                self.deg = -1       # when we`ll add an element to list of koeffs deg`ll be increased


        def eq(self, G):
            if G.deg != self.deg:
                return False
            for i in range(0, self.deg + 1):
                if self.koeffs[i].f_x != G.koeffs[i].f_x:
                    return False
            return True 

        # updates degree and also deletes heading 0*y^t:  0*y^3 + x*y^2 + 1 ---> x*y^2 + 1  aka [1, x, 0] ---> [1, x]   
        def update_deg(self):
            deg = len(self.koeffs) - 1
            if deg == -1: print("!!! Invalid polynom(not inited)"); self.deg = -1                      # if poly is undefinded(no koeffs) => return invalid value

            for koef_x in self.koeffs[::-1]:     # try to find first y^deg with not zero koefficent by it
                if koef_x.f_x != 0:
                    self.deg = deg
                    break
                self.koeffs.pop()               # delete empty(0) koefs by y of highest degree. Ex: otherwise 0*y^3 + y^2 - unintuitive => delete 0*y^3
                deg -= 1
            else:
                self.deg = 0                        # all koeffs = 0 => self = 0
                self.koeffs = [poly_Z_2(f_x=0)] # all koeffs were poped => koefs = [], but must be [0] - we`ve added here

        def __mul__(self, G):
        
            mul = poly_Y(GF=self.GF)
            for deg in range(0, self.deg + G.deg +1):   
                mul.koeffs.append(poly_Z_2())                   # giving memory for mul[self.deg, self.deg + G.deg]

            for deg in range(0, self.deg + 1):                  # +1 because ammount of koeffs = degree of poly + 1 
                if self.koeffs[deg].f_x != 0:  
                    for G_deg in range(0, G.deg + 1):
                        if G.koeffs[G_deg].f_x != 0:
                            #print(f"    +++ to deg = {deg + G_deg} | {self.koeffs[deg]} * {G.koeffs[G_deg]} = {self.GF.mul(self.koeffs[deg], G.koeffs[G_deg])}")
                            mul.koeffs[deg + G_deg] += self.GF.mul(self.koeffs[deg], G.koeffs[G_deg])     # cheap multiplying from GF`s table

            mul.update_deg()      
            return mul



        def __add__(self, G):
            sum = poly_Y(self.GF)             
            
            if self.deg >= G.deg:
                for i in range(0, G.deg + 1):               # self and G poly may be of different degree - so summ `em while they have same degrees
                    sum.koeffs.append((self.koeffs[i] + G.koeffs[i]) )
                for i in range(G.deg + 1, self.deg + 1):
                    sum.koeffs.append(self.koeffs[i]) 

            else: 
                for i in range(0, self.deg +1):            # self and G poly may be of different degree - so summ `em while they have same degrees
                    sum.koeffs.append((self.koeffs[i] + G.koeffs[i]) )
                for i in range(self.deg +1, G.deg +1):
                    sum.koeffs.append(G.koeffs[i]) 

            sum.update_deg()                            # modulo could decreaze power from self.deg + G.deg to someth lower
            return sum


        # Итерратор по poly_Y
        def next(self):
            F_y = deepcopy(self)

            for i in range(0, len(F_y.koeffs)):
                if F_y.koeffs[i].f_x != 2**(self.GF.p_x.deg) - 1:     # check if koef isn`t max (in GF(2^t)) and can be increased            
                    F_y.koeffs[i] = F_y.koeffs[i].next()
                    for j in range(0, i):
                        F_y.koeffs[j] = poly_Z_2()
                    return F_y

            else:                                   # if all koeffs are max
                for i in range(0, len(F_y.koeffs)):
                    F_y.koeffs[i] = poly_Z_2()           # zeroing all the koeffs
                F_y.koeffs.append(poly_Z_2(f_x = 1))     # append new(and minimal) koeff to the end
                F_y.deg += 1
                return F_y

        # return y^k
        def y_k(self, k: int):
            y_k = poly_Y(GF=self.GF)    
            for deg in range(0, k):
                y_k.koeffs.append(poly_Z_2())   # append k - count zero koeffs
            
            y_k.koeffs.append(poly_Z_2(f_x=1))  # one more element that`s = 1
            y_k.deg = k
            return y_k



        # For poly_Y: return modulo(self % G) and devider(self/G) 
        def mod_dev(self, G):
            
            def sub(res:poly_Y, div: poly_Y):
                k = res.deg - div.deg
                dev = self.y_k(k)                                                   # y^k
                dev.koeffs[-1] = self.GF.mul(res.koeffs[-1],  self.GF.get_opposite(div.koeffs[-1]) )   # biggest koef by res and div_i must be equal to continue deviding => dev = res[-1]/div[-1]. ps /div[-1] <=> * div[-1]^-1 
                
                div_i = div * dev                                                   # and here dev[-1] * div[-1] <=> div_i[-1] = res[-1]  
                
                return res + div_i, dev             # res - div  is same as  res + div in GF(2^t)
                    

            res = deepcopy(self)
            dev = poly_Y(self.GF, koeffs=[poly_Z_2(f_x = 0)])   # 0 * y^0 - 0 as F(y)

            while res.deg > G.deg:                
                res, dev_i = sub(res, G)
                dev += dev_i
                
            # Ex: res = 0, dev_i = alpha^9 - this might cycle forewer above(due no changes and condition remains true) - so I separated this case 
            if res.deg == G.deg:   
                res, dev_i = sub(res, G)
                dev += dev_i

            return res, dev 


        def on_val(self, F: poly_Z_2):
            Val = poly_Z_2()
            deg = 0
            for f in self.koeffs:
                if f.f_x != 0:
                    Val += self.GF.mul(f, self.GF.pow(F, deg))
                deg += 1
            return Val



        def __str__(self):

            f_y_str = ""
            k_deg = 0
            for koef in self.koeffs:
                if koef.f_x != 0:  f_y_str += f" α^{self.GF.get_ord(koef)}*y^{k_deg}  + "     # ({koef}) alph^{self.GF.get_ord(koef)}
                k_deg += 1

            # pretty output for 1-th and 0-th powers
            f_y_str = f_y_str.replace("y^1 ", "y")
            f_y_str = f_y_str.replace("*y^0 ", "")

            return f_y_str.rstrip(" + ") 


# Some exceptions:
class Missfit_t_and_m(Exception):
    def __init__(self, t, m):        
        # переопределяется конструктор встроенного класса `Exception()`
        sep = " \n" + "!"*40 + "\n"
        super().__init__(f"{sep} 2^{m} - 1 < minimal code distance({2*t + 1}) - u should increase m too fix t={t} errors{sep}")

class Too_small_m(Exception):
    def __init__(self, m):
        sep = "\n" + "!"*40 + "\n"
        super().__init__(f"\n m = 1 is too small{sep}")

class Error_too_big_message(Exception):
    def __init__(self, Message, INFO_BLOCK):
        sep = "\n" + "!"*40 + "\n"
        super().__init__(f"\n{sep} len(Message) in bits = {len(Message)}*8 = {len(Message)*8},\nwhile INFO_BLOCK size = {INFO_BLOCK} bits {sep}")

class BAD_Infoblock(Exception):
    def __init__(self, INFO_BLOCK):
        sep = "\n" + "!"*40 + "\n"
        super().__init__(f"\n{sep} INFO_BLOCK is too small = {INFO_BLOCK} bits {sep} for chosen m. Please, increase m")



class BCH:
    
    # t - ammount of errors fixing,  m - defines size of block = 2**m -1 = |GF(2**m)/{0}|
    def __init__(self, t: int, m: int):    
        self.t = t

        if m < 2:
            raise Too_small_m(m)

        # minimal code distance >= 2t + 1  --->  n=2^m-1  >= 2*t + 1 
        if 2**m - 1 < 2*t + 1:
            raise Missfit_t_and_m(t,m)
        
        # if deg(p(x))=m <= 10 works fast 
        p, GF_tbl = poly_Z_2.get_primitive(m)       # find poly that can permutate all GF(2^m)/{0} elements. 
        self.GF_x = GF2_x(p, Table=GF_tbl)          # build GF(2^m)

        # g`s roots = GF_x.table[0:delta] : alpha^0 ... alpha^(2t+1) - roots of g(x). delta = 2t+1 - minimal code distance 
        self.g_x = self.GF_x.get_g_x(t)             # g(x) - generator. p(x) - generator. because it`s primitive poly(permutes every != 0 poly of deg < m)
        
        self.CODED_BLOCK = 2**m - 1 
        self.INFO_BLOCK = self.CODED_BLOCK - self.g_x.deg
        if self.INFO_BLOCK == 1 and m > 3:
            raise BAD_Infoblock(self.INFO_BLOCK)

        sep = "\n" + "=" * 50 + "\n"
        print(f"{sep}  Built BCH-coder that fixes t = {t} errors \n   INFO_BLOCK = {self.INFO_BLOCK} bits, \n   CODED_BLOCK = {self.CODED_BLOCK} bits\n   Min. code distance >= {2*t + 1}\n   Safety = {t/self.CODED_BLOCK} (t/CODED_BLOCK)\n   Rate = {self.INFO_BLOCK/self.CODED_BLOCK} (inf/code block) {sep}") 

    

    # some chanel distortion
    def distortion(self, k: int):

        dirty = poly_Z_2()        

        for _ in range(0, k): 
            i = random.randint(0, len( bin(self.CODED_BLOCK).lstrip("0b")) )
            dirty.f_x += 2**i

        dirty.update_deg()

        return dirty


    # messege --> poly-blocks 
    def str_encode(self, message: str):
        mess = bytes(message, encoding = "ASCII")
        enc = ""

        # bytes --> binary string
        for ch in mess:
            #print(f" {ch} {bin(ch)}")
            enc += '{:0>8b}'.format(ch)         # Adding 8 charectered string with heading(>) zeroes(0) for binary(b) representation of ch.  see "format_spec" - syntax in format()

        print(f"Message in binary: \n{enc}\n\n{chr(9935)} Now devide it into blocks and convert`em into polynoms:\n")

        # slicing binary string into blocks of INFOBLOCK size
        blocks_x = []
        i = 0

        if len(enc) < self.INFO_BLOCK:
            blocks_x.append( poly_Z_2(str_koeffs=enc) )
            return  blocks_x


        while i < len(enc):
            if len(enc) - i - self.INFO_BLOCK >= 0:
                print(f"+{enc[i:i+self.INFO_BLOCK]} | {poly_Z_2(str_koeffs=enc[i:i+self.INFO_BLOCK])} ")
                blocks_x.append( poly_Z_2(str_koeffs=enc[i:i+self.INFO_BLOCK]) )    
            else:                                                                             # in case last block can`t be of size INFO_BLOCK 
                print(f"-{enc[i:i + len(enc) - self.INFO_BLOCK]} | {poly_Z_2(str_koeffs=enc[i: i + len(enc) - self.INFO_BLOCK])} ")
                blocks_x.append( poly_Z_2(str_koeffs=enc[i:i + len(enc) - self.INFO_BLOCK]))     # len(enc) - self.INFO_BLOCK - ammount of remained bits 
            i += self.INFO_BLOCK
 
        return blocks_x


    # poly-blocks ---> messege
    def str_decode(self, blocks_x:list):
        
        bin_repr = ""
        for F in blocks_x:                            
            bin_repr += bin(F.f_x).lstrip("0b")[::-1]                                   # integers in polynoms stores in reversed order 
            bin_repr += "0" * (self.INFO_BLOCK - len(bin(F.f_x).lstrip("0b")[::-1]))    # padding to length of info block each block(ex: 2 = 10 ---> 0010 if inf_block = 4)

        bin_repr = bin_repr + "0" * (len(bin_repr) % 8)   # padding to byte sized str. ps byte = 8 bits
        
        dec_str = ""
        i = 0
        #print("V: " + bin(ord("V")))
        while i < len(bin_repr):
            #print("  " + bin_repr[i: i+8])
            dec_str += chr( int("0b" + bin_repr[i: i+8], base=2))
            i += 8
        

        return dec_str



    # syndrom-polynom(y): each i-th koeff = value of recieved code-vector(v`) on i-th power of alpha from GF_X
    def get_S_y(self, recived: poly_Z_2):       # recived - v`(x) - code-message with some distortion from chanel 

        S_y = poly_Y(self.GF_x)                 
        koefs = recived.to_list()               # getting representation of polynom as list of it`s koeffs
                
        # give memory for each S(y)`s koeff
        for deg in range(0, 2 * self.t):
            S_y.koeffs.append(poly_Z_2())

        deg = 0                                 # degree of y
        for k in koefs:
            if k != 0:
                for i in range(0, 2 * self.t):
                    S_y.koeffs[i] += self.GF_x.pow(self.GF_x.get_alpha_i(deg), i+1)   # Updating all koeffs of S(y): S(y)_i += (alpha^deg)^(i+1). ps i+1 bc i migh be 0   
            deg +=1
        
        
        S_y.update_deg()
        
        '''
        print("Let`s find S(c`(x)):")
        for S_i in S_y.koeffs:
            if S_i.f_x != 0:
                print(f"S_i = {S_i} <=> α^{self.GF_x.get_ord(S_i)}")
        print(f"\nS(y) = {S_y}\n")
        '''

        return S_y                              # S(y) = alpha^k_1 + alpha^k_2 * y^1 + ... + alpha^k_2t(x) * y^2t-1  - this poly(y) stores as list of koeffs(has type f(x))     


    
    # apaptation(reducing some steps) of extended Euclid`s algorithm for BCH`s purposes(finding sigma(y))
    # deg a(x) >= deg b(x)
    def reduced_ext_GCD(self, a: poly_Y, b: poly_Y):
        
        A_y, B_y  = deepcopy(a), deepcopy(b)                   # create physical coppies of givven polynoms not to corrupt `em`
        q, r      = poly_Y(self.GF_x), poly_Y(self.GF_x)
        c_b, c_b_ = poly_Y(self.GF_x, koeffs=[poly_Z_2(f_x=0)]), poly_Y(self.GF_x, koeffs=[poly_Z_2(f_x=1)])      # 0 and 1 as polynoms
        
        zero_Y = poly_Y(self.GF_x, koeffs=[poly_Z_2()])                 # 0  in poly_y type


        while (not B_y.eq(zero_Y)) and A_y.deg >= self.t:
            print(f"A_y: {A_y}, B_y: {B_y}")
            r, q = A_y.mod_dev(B_y)
            print(f"    q: {q}  r: {r}")

            A_y,B_y = B_y, r 
            c_b, c_b_ = c_b_, c_b + q * c_b_   

            print(f"    c_b: {c_b}\n")

        return c_b      # almost sigma(y). koeff by sigma by def is 1, so we can multiply on inverted koeff[0] all the sigma(y), but it`s roots`ll remain same, so it`s useless

    # messege as poly-blocks ---> BCH-encoded poly-blocks
    def encode(self, blocks_x: list):
        print("\nCoding each poly-block to BCH poly:")
        enc_blocks = []
        for F in blocks_x:
            print(f"  {F*self.g_x}")
            enc_blocks.append(F*self.g_x)

        return enc_blocks

    def is_ok(self, koeffs, val):
        for el in koeffs:
            if el.f_x != val:
                return False
        return True


    # BCH-decoder - for single block. Multiple blocks - use Terminal
    # BCH-encoded poly ---> messege poly 
    def decode(self, v_: poly_Z_2):

        e = poly_Z_2()                  # polynom of errors 
        S_y = self.get_S_y(v_)          # syndrom(y)

        if self.is_ok(S_y.koeffs, 0):   # if no errors
            print("There`s no errors")
            return deepcopy(v_/self.g_x)                        


        print(f"\nLet`s find sigma(y) as red_ext_GCD(S(y), y^k):\n   S(y) = {S_y}\n   y^k = {poly_Y(self.GF_x).y_k(2 * self.t)}\n")

        sigma_y = self.reduced_ext_GCD(poly_Y(self.GF_x).y_k(2 * self.t), S_y)  # red_GCD(y^2y, S_y) => locator 

        print(f"sigma(y) = {sigma_y}")

        #print(f"\nS(y) = {S_y}\ny^k = {poly_Y(self.GF_x).y_k(2 * self.t)}\nsigma(y) = {sigma_y}")

        # locator aka sigma(y) by def: sigma(y) = (1 + a^j1 * y)(1 + a^j2 * y)...(1 + a^jr * y), where j - is position of error
        err_counter = 0
        for alpha in self.GF_x.table:
            if sigma_y.on_val(alpha).f_x == 0:                      # locator(alpha) = 0 <=> alpha is root =>  alpha^-1 is in e(x) 
                
                err_counter += 1                                    # new error find
                opos = self.GF_x.get_opposite(alpha)                # alpha^-1 = alpha^j => j = get_ord(alpha^-1) 
                e += poly_Z_2(f_x = 2 ** self.GF_x.get_ord(opos))   # v += x^j
                print(f"    sigma(α^{self.GF_x.get_ord(alpha)}) = 0 => {self.GF_x.get_ord( self.GF_x.get_opposite(alpha) )}  bit is corupted")


        if err_counter > self.t:
            print(f"\nThere`s more than {self.t}(max) errors \\_(0_0)_/ ")
            return 0
        else:
            print(f"\nFound {err_counter} error(s), e(x) = {e}")
            return (v_ + e) /self.g_x                              # reduce e from v_ .ps + and - same GF(2)
    

def Terminal(message: str, t_errors: int, m: int):
    print(f"Sending message: {message}")

    Coder = BCH(t_errors, m)
    print(f"\nGenerator:\ng(x) = {Coder.g_x}\n")
    v_blocks = Coder.str_encode(message)
    enc_blocks = Coder.encode(v_blocks)

    k = int(input(f"{chr(10067)} How many errors happens while message transition(in each block)???\n : "))
    lst_e_x = [Coder.distortion(k) for block in enc_blocks]                                # forming distortion for each block

    dec_blocks = []
    for i in range(0, len(enc_blocks)):
        print(f"\n=== Decoding {i}-th block:")
        dec_blocks.append( Coder.decode(enc_blocks[i] + lst_e_x[i]) )

    print(f"\nLet`s decode(coder fixes {t_errors} - max), g(x) = {Coder.g_x}")


    dec_str = Coder.str_decode(dec_blocks)

    magic = "      " + chr(1161)*30
    print(f"\n{magic*7}\n       Decoded message:\n        {dec_str}\n{magic*7}")


    ans = input(f"{chr(10067)} Wanna see used GF?[y]\n:")
    if ans == "y":
        print(Coder.GF_x)



if __name__ == "__main__":
    
    
    GF2_5 = GF2_x(poly_Z_2(koeffs=[1, 0, 1, 0, 0, 1]))       # x^5 + x^2 + 1 - primitive poly 
    GF2_4 = GF2_x(poly_Z_2(koeffs=[1, 1, 0, 0, 1]))       # x^4 + x + 1 - primitive poly


    t = int(input("Enter t - ammount of errors coder must to decode\n :"))
    m = int(input("Enter m:  len(coded block) = 2^m - 1 = |GF(2^m)\\0|\n :"))
    message = input("Input your message\n :")
    Terminal(message, t, m)



    '''
    #F_y = poly_Y(GF2_4, koeffs = [poly_Z_2(f_x=3), poly_Z_2(f_x=7), poly_Z_2(f_x=2) ])
    #F_y_2 = poly_Y(GF2_4, koeffs = [poly_Z_2(f_x=2), poly_Z_2(f_x=6)])
    #F_y_3 = poly_Y(GF2_4, koeffs=[poly_Z_2()])

    Coder = BCH(2, 4)
    print(f"g(x) = {Coder.g_x}")
    v = poly_Z_2(str_koeffs="1010101")
    coded_v = v * Coder.g_x
    e_x = poly_Z_2(str_koeffs="10001")

    print(f"\nv(x) = {v} ---> c(x) = {coded_v}\ne(x) = {e_x}\n\nc`(x) = {coded_v + e_x} \n- Decoder must find it\n")

    v_ = Coder.decode(coded_v + e_x)
    print(f"c`(x) ---> {v_} - decoded")
    '''

    '''
    GF = GF2_4
    
    A = poly_Y(GF).y_k(6)                                # y^6
    B = poly_Y(GF, koeffs_alpha=[7, 14, 13])  
    
    coder_ = BCH(2, 4)
    print(f"{coder_.GF_x}")
    sigma_y = coder_.reduced_ext_GCD(A, B)
    

    # print(f" GCD [({A}),({B})]:\n   sigma(y): {sigma_y}\n")
    '''


    '''
    Coder = BCH(7, 5)
    F = Coder.str_encode("012345678910111213141516171819202122232425262728293031323334353637383940012345678910111213141516171819202122232425262728293031323334353637383940012345678910111213141516171819202122232425262728293031323334353637383940012345678910111213141516171819202122232425262728293031323334353637383940012345678910111213141516171819202122232425262728293031323334353637383940012345678910111213141516171819202122232425262728293031323334353637383940")
    Coded = Coder.encode(F)

    F_dec = []
    for i in range(0, len(Coded)):
        F_dec.append(Coder.decode(Coded[i] + Coder.distortion(7)))
    
    for i in range(0, len(F_dec)):
        if F[i].f_x != F_dec[i].f_x:
            print(f"i = {i}| Decoded: {F_dec[i]}\nOriginal: {F[i]}")

    print(f"\nF: {len(F)}, F_dec: {len(F_dec)}\nDecoded: {Coder.str_decode(F_dec)}")
    '''