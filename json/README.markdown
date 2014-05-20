README
======

Matlab JSON library

This package contains Matlab class to serialize/decode matlab object in
json format. The software uses org.json java package to convert json to
java object and then translates it into Matlab object.

Kota Yamaguchi 2012 <kyamagu@cs.stonybrook.edu>


Usage
-----

Make sure your matlab path includes a directory containing @JSON. Matlab
recognizes a directory whose name starts with '@' as a class directory.

To serialize matlab object:

    >> X = struct('matrix', magic(2), 'char', 'hello');
    >> S = JSON.dump(X);
    >> disp(S);
    {"char":"hello","matrix":[[1,3],[4,2]]}

To decode json string:

    >> X = JSON.load(S);
    >> disp(X);
          char: 'hello'
        matrix: [2x2 double]


Limitation
----------

 * Currently the software doesn't support N-D array.
 * Due to the multiple ways to represent an array in Matlab (i.e., numeric
   array, cell array, or struct array), it is not possible to represent
   everything in compatible format. The software converts cell array of
   numbers or struct to numeric array or struct array when possible.

License
-------

You may redistribute this software under BSD license.


See also
--------

JSON in Java: http://json.org/java/