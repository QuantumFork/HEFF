︠c10b506e-4bfd-449e-80b8-3cf3311a8ebcs︠
def divisor_eea(a,b,N): #When xgcd(a,b) fails, it fails because you have an element which is not
    #invertible mod N at some step. Taking the gcd of that element with N gives a non-trivial
    # divisor, which is how HECF factorization works.
    storea = a
    storeb = b
    P.<x> = Integers(N)[]
    l = int(P(b).leading_coefficient())
    g = gcd(l,N)
    if g != 1 and g!= N: #This is one place where eea might fail, outputting a non-trivial divisor
        return g, b,0
    else:
        while b!= 0:
            l =int(P(b).leading_coefficient())
            g = gcd(l,N)
            if g != 1 and g!= N:
                return g, b,0 
        # By going through, what is essentially the euclidean algorithm, and checking if 
        # Each step is possible, we ensure that xgcd(storea,storeb) will work.
        # If it won't it will be because our code has found a non-invertible element of Z/NZ.
            else:
                r=a%b
                a=b
                b=r
        return xgcd(storea,storeb) #By this point we are sure xgcd() will work, so we can safely use it.



def reduce_divisor(a,b, f, h, N):
    # The output of the operations of addition and multiplication will not, in general, be semi-reduced divisors.
    #This code reduces them. 
    C = HyperellipticCurve(f,h)
    genus = C.genus()
    while P(a).degree() > genus: #The code follows Cantor's algorithm exactly.
        reduceda = (f -b*h-b^2)//a
        lca = P(reduceda).leading_coefficient()
        g1 = gcd(lca,N) 
        if g1 != 1 and g1 !=N:# Since we are not, in general, working over a field, g1 may be a factor of N
                              # If it is, the code outputs a non-trivial divisor of N.
            return (g1, 0)
        reducedb = (-h-b)%reduceda
        a = reduceda
        b = reducedb
        a = a/P(a).leading_coefficient()
    return a,b




def add_divisors(a1,b1,a2,b2, f, h,N): #point addition in a hyperelliptic curve
    if gcd(divisor_eea(a1,a2,N)[0],N) !=1 and P(divisor_eea(a1,a2,N)[1]).degree() == 0: 
        #this simply checks if divisor_eea output a non-trivial divisor of N.
        return divisor_eea(a1,a2,N)[0], 0
    else: #The rest is just cantors algorithm.
        d1, e1, e2 = divisor_eea(a1,a2,N) 
        if d1 == 1: # This if statement is part of Cantor's algorithm.
            a = a1*a2
            b = (e1*a1*b2 + e2*a2*b1)%a
        else:
            d, c1, s3 = divisor_eea(d1, b1 + b2 + h,N)
            s1 = c1*e1
            s2 = c1*e2
            if P(d^2).degree() == 0 and gcd(d^2, N) !=1: 
            #This is another place where the algorithm might fail, and we can extract a non-trivial divisor.
                return gcd(d^2, N), 0
            a=a1*a2//d^2
            if P(d).degree() == 0 and gcd(d, N) != 1 and gcd(d,N) !=N: 
            #This is another place where the algorithm might fail, and we can extract a non-trivial divisor.
                return gcd(d,N), 0
            b=((s1*a1*b2+s2*a2*b1+s3*(b1*b2+f))//d)%a
        return reduce_divisor(a,b,f,h,N)



def multiply_divisor(a,b, k, f, h, N):
    #point multiplication in a hyperelliptic curve, which uses the double-and-add algorithm.
    sumu = a
    sumv = b
    binary_k = list(format(k, "b"))
    l = len(binary_k)
    for i in range(1,l):
        if int(binary_k[i]) == 0:
            if P(sumu).degree()==0:
                if gcd(sumu, N) !=1:
                    return sumu, sumv
            sumu, sumv = add_divisors(sumu, sumv, sumu, sumv, f, h, N)
        if int(binary_k[i]) == 1:
            if P(sumu).degree()==0:
                if gcd(sumu, N) !=1:
                    return sumu, sumv
            sumu, sumv = add_divisors(sumu, sumv, sumu, sumv, f, h, N)
            if P(sumu).degree()==0:
                if gcd(sumu, N) !=1:
                    return sumu, sumv
            sumu, sumv = add_divisors(sumu, sumv, a, b, f, h, N)
    return sumu, sumv


def get_point(N, f):
    #This allows us to generate a point on an elliptic curve of genus 2, with h = 0, to be used hef().
    y = polygen(Zmod(N))
    R = y.parent()
    for i in (0..(N-1)//2):
        r = Mod(f(i),N)
        if r.is_square():
            return(y-i), R(sqrt(r))

def compute_k(B): 
    #The hyperelliptic curve factorization algorithm uses a bound k. These are good numbers
    #To multiply a divisor by in order to try to implement hecf.
    k = 1
    p = Primes()[0]
    while p <= B:
        n = 1
        while p^(n+1) <= B:
            n+=1
        k*= p^n
        p = Primes().next(p)
    return(k)


def hecf(N, B, f):
    #This implements 1 iteration of HECF. If this outputs a point, HECF requires you to try a new curve.
    #I often found it's easier to just choose a larger B.
    p.<x> = Integers(N)[]
    h = P(0)
    f = x^5 + 3*x + 40 #This choice of hyperelliptic curve worked for all inputs I tried.
    #In later versions of this code, one could remove this, and allow for an arbitrary f
    # However, it's not necessary here.
    a1,b1 = get_point(N,f)
    k = compute_k(B)
    m = multiply_divisor(a1,b1,k,f,h,N)
    if P(m[0]).degree()==0:
        if gcd(m[0],N) !=1:
            return "Computation of kP, failed. A non-trivial factor of", N, "is", m[0]
    else:
        return k, "P =", m

N = 1234567891011
B = 2
P.<x> = Integers(N)[]
f = x^3 + 5*x + 40
h = P(0)

hecf(N,B, f)
︡0df5c57d-526c-48d2-9530-234e2c1dec60︡{"stdout":"('Computation of kP, failed. A non-trivial factor of', 1234567891011, 'is', 13)\n"}︡{"done":true}









