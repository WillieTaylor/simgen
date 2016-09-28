#!/usr/bin/env python
### hilbert_test.py -- tests for hilbert.py


from string import join
from sys import argv
from math import log, ceil
from copy import copy

from hilbert import gray_encode_travel, gray_decode_travel, child_start_end
from hilbert import int_to_Hilbert, Hilbert_to_int



def n2digits( n, modulus, nd ):
    return [ ( n / modulus**p ) % modulus for p in range( nd-1, -1, -1 ) ]


def bin2str( n, nbits ):
    return "".join( [ chr( ord('0') + digit ) for digit in n2digits( n, 2, nbits ) ] )


def str2bin( str ):
    return reduce( lambda n,bit: 2*n+bit, [ ord(c)-ord('0') for c in str ] )


def gray_test():
    for n in range( 16 ):
       g = gray_encode(n)
       assert gray_decode(g) == n

       print bin2str( n, 4 ), bin2str( g, 4 )



def hilbert_encode_partway_test( nbits, ndims ):
    for n in range( 2 ** ( nbits * ndims ) ):
        encoded = hilbert_encode_partway( n, nbits, ndims )
        for digit in encoded:
            print bin2str( digit, nbits ),
        print




def gray_code_travel_test( nbits ):
    encode = gray_encode_travel
    decode = gray_decode_travel
    mask = 2 ** nbits - 1
    onebit = [ 1 << i for i in range( nbits ) ]
    for i in range( 2 ** nbits ):
        for start in [ 0, mask ]:
            for sh in range( nbits ):
                offs = 1 << sh
                end = start ^ offs
                g = encode( start, end, mask, i )
                if i == 0:
                    assert g == start
                else:
                    prev_g = encode( start, end, mask, i-1 )
                    assert ( prev_g ^ g ) in onebit  # One bit changed.
                if i == mask:
                    assert g == end
                j = decode( start, end, mask, g )
                assert j == i
                print bin2str( g, nbits ),
            print "  ",
        print


def stars( n, nbits ):
    return "".join( [ " *"[digit] for digit in n2digits( n, 2, nbits ) ] )


    


def child_start_end_test( nbits ):
   mask = 2 ** nbits - 1
   parent_start = 0
   parent_end = 2 ** ( nbits - 1 )  # Canonical Gray code end.

   heads = [ "i", "parent", "child" ] 
   format = "%%%ds  %%%ds  %%%ds" % tuple( [ max(nbits,len(s)) for s in heads ] )
   print format % tuple( heads )

   for i in range( mask + 1 ):
       parent = gray_encode_travel( parent_start, parent_end, mask, i )
       child_start, child_end = child_start_end( parent_start, parent_end, mask, i )
       child = child_start
       if i > 0:
           parent_xor = parent ^ parent_prev
           child_xor = child ^ child_prev
           print format % ( " ", stars( parent_xor, nbits ), stars( child_xor, nbits ) ),
           if parent_xor != child_xor:
               print "    BZZT!",
           print
       print format % ( "%d"%i, bin2str( parent, nbits ), bin2str( child, nbits ) ),
       if i == 0:
           if child != parent:
               print "    CHILD != PARENT !!",
       print
       child = child_end
       print format % ( " ", " ", bin2str( child, nbits ) ),
       if i == mask:
           if child != parent:
               print "    CHILD != PARENT !!",
       print
       parent_prev, child_prev = parent, child



def child_start_end_tests( ):
    for nbits in range( 1, 5 ):
        print
        print
        child_start_end_test( nbits )


## check_visited -- Check that every point in the cube was visited.
#
# visited is a dictionary (or "hash" or associative array) whose key
# is a tuple representing the point coordinates.  Prefix is a list
# of coordinates that's built up by recursion to the right number of
# dimensions.  "tuple(prefix) in visited" converts the list to a tuple 
# and asks whether an entry with that key is in the dictionary.
def check_visited( visited, prefix, linear_size, nD ):
    if nD > 0:
        for i in range( linear_size ):
            check_visited( visited, prefix + [i], linear_size, nD - 1 )
    else:
        if not tuple( prefix ) in visited:
            print "missing", prefix


def tuple_print_width( n ):  # tuple with n 1-digit numbers and trailing sp
    if n == 0:  return 3   # ()_
    if n == 1:  return 5   # (1,)_  special case
    else:       return n * 3 + 1


def int_to_Hilbert_test( n, nD, doPrint=True ):
    tups_per_line = 2**( int( log( 80 / tuple_print_width( nD ), 2 ) ) )
    linear_size = int( ( n + .5 ) ** (1./nD) )
    isCube = linear_size ** nD == n
    if not isCube:
        print "(This is not a complete cube:)"
    visited = {}        
    for i in range( n ):
        pt = tuple( int_to_Hilbert( i, nD ) )
        if doPrint: print pt,
        visited[ pt ] = True      # Add pt to "visited" dictionary.
        j = Hilbert_to_int( pt )
        if j != i:
            print "BZZT, decodes to", j,

        if i > 0:  # Check that we have stepped one unit, along one dim:
            nDiffs = sum( [ pt[k] != prev_pt[k] for k in range( nD ) ] )
            if nDiffs != 1:
                print "BZZT, different along", nDiffs, "dimensions",
            maxDiff = max( [ abs(pt[k] - prev_pt[k]) for k in range( nD ) ] )
            if maxDiff != 1:
                print "BZZT, max difference is", maxDiff,

        prev_pt = pt
        if ( i + 1 ) % tups_per_line == 0:
            if doPrint: print

    if isCube:
        check_visited( visited, [], linear_size, nD )


if argv[0].endswith( "/hilbert_test.py" ) or argv[0] == "hilbert.py":
    # gray_code_travel_test( 4 )
    # child_start_end_tests()
    int_to_Hilbert_test( 8, 1 )
    print
    int_to_Hilbert_test( 64, 2 )
    print
    int_to_Hilbert_test( 64, 3 )
    print
    int_to_Hilbert_test( 4096, 3, doPrint=False ) # nChunks > nD
    print
    int_to_Hilbert_test( 4096, 4, doPrint=False )
    print
    k = 10**12  # definitely a long
    for nD in [ 2, 3, 4, 5 ]:
        # Decode point ( k, 0, 0... ) to an int...how many digits?:
        pt = [ 0 ] * nD
        pt[0] = k
        x = Hilbert_to_int( pt )
        print x, log( x, 10 )
        # Double-check:
        pt2 = list( int_to_Hilbert( x, nD ) )
        if pt2 != pt:
            print "BZZT!  Encodes to ", pt2

        # Encode k to an nD point:
        pt = int_to_Hilbert( k, nD )
        print pt
        # Double-check:
        y = Hilbert_to_int( pt )
        if y != k:
            print "BZZT!  Decodes to", y


