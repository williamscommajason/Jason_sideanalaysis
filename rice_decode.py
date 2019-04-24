import struct
import os

def BitString(BitString):

    for bit in BitString:
        yield bit


def decompress(f):

    
    #with open(filename, "rb") as f:

    lists = []
    byte = ""
    while f.getbuffer().nbytes > f.tell(): 

        k = struct.unpack('@i',f.read(struct.calcsize('i')))[0]
        size = struct.unpack('@Q',f.read(struct.calcsize('Q')))[0]         
        bString = ""

        for i in range(size): 

            byte = f.read(1)
            hexbyte = struct.unpack("@B",byte)[0]
            binary = "{0:b}".format(hexbyte)

            if len(binary) < 8:
                for i in range(8 - len(binary)):
                    binary = '0' + binary

            bString = bString + binary
                             
    
        codes = decode_bitString(bString,k)
        rice_dictionary = rice_dict(k,50)
   
        unsigned = []
    

        for i in codes:
            if i in rice_dictionary.keys():
                unsigned.append(rice_dictionary[i])
            else:
                unsigned.append(decode_rice_byte(i,k))

        signed = back_to_signed(unsigned)

        lists.append(signed)    
             
    return lists
            

def decode_bitString(bString, k):
 
    zero_flag = False
    one_flag = False
    codes = []
    
    bString = BitString(bString)

    while True:

        bit = next(bString)

        zero_flag = False
        one_flag = False
        quotient = ''
        remainder = ''
        try:
            
            if bit == '0' and one_flag == False:
                
                quotient = '0' 
                zero_flag = True

            if zero_flag == True and one_flag == False:
                for i in range(k-1,-1,-1):

                    bit = next(bString)
                    remainder = remainder + bit    

                codes.append(quotient + remainder)
                zero_flag = False
                continue
        
            while bit == '1' and zero_flag == False:

                one_flag = True
                quotient = quotient + '1'
                bit = next(bString)

                if bit == '0':

                    zero_flag = True
                    quotient = quotient + bit
                    break
        
               
            if zero_flag == True and one_flag == True:              
                for i in range(k-1,-1,-1):

                    bit = next(bString)
                    remainder = remainder + bit    
                codes.append(quotient + remainder)

        except StopIteration:
            #print("Last element was: ", bit)
            break
        
    return codes
            
def decode_rice_byte(rb,k):
    
    rb = BitString(rb) 
    zero_flag = False
    one_flag = False
    
    while True:
    
        try:
            bit = next(rb)
            quotient = ''
            if bit == '0' and one_flag == False:
                for i in range(k-1,-1,-1):
                    bit = next(rb)
                    quotient = quotient + bit

                x = int(quotient,2)
                
                return x

            ones = 0   
            while bit == '1':
                one_flag = True
                ones += 1  
                bit = next(rb)
                if bit == '0':
                    zero_flag = True
                    quotient = quotient + bit
                    break
 
            if zero_flag == True and one_flag == True:
                remainder = ''
                for i in range(k-1,-1,-1):
                    bit = next(rb)
                    remainder = remainder + bit
                    

            x = ones*(2**k) + int(remainder,2)
           
            return x

        except StopIteration:
            break  
    
    
def rice_dict(k,end):

    rice_dictionary = dict()
    ints = tuple(range(end+1))

    for x in ints:
        q = x / (1 << k)
        q = int(q)
        bString = ''
        for i in range(q):
            bString = bString + '1'
        bString = bString + '0'
        for i in range(k-1, -1, -1):
            bString = bString + str((x >> i) & 1)
            
        rice_dictionary[bString] = x

    return rice_dictionary
  
def undiff(L):

    for i in range(1,len(L)):
        L[i] += L[i-1]
    return L
                   
def back_to_signed(L):

    signed = []
    for x in L:
        if x%2 == 0:
            signed.append(int(x/2))
        if x%2 == 1:
            signed.append(-int((x+1)/2))

    return signed
        

if __name__ == "__main__":

    f = open('rice.bin', 'w+b')
    signed = decompress(f)
    print(rice_dict(2,20))
    
    print(signed)
    
