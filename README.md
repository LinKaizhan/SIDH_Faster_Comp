This code is complementary to the IACR eprint paper .

we mainly improve the performance of the public-key compression of SIDH, especially the efficiency and the storage of pairing 
computation involved. Our experimental results (Method1) show that the memory requirement for pairing computation is reduced
by a factor of about 1.31, and meanwhile, the instantiation of the key generation of SIDH is 3.99%~ 5.95% faster than the current 
state-of-the-art. Besides, in the case of Bob, we present another method (Method 2) to further reduce storage cost, while the 
acceleration is not as obvious as the former. 

PQCrypto-SIDH-master_M1 is the code when Bob utilizes Method 1 for pairing evaluation,  and PQCrypto-SIDH-master_M2 is the 
code when utilizing Method 2.

Our code is based on SIDH v3.3 (https://github.com/Microsoft/PQCrypto-SIDH), and the usage of our code is as same as the latter.
See README.md in the file PQCrypto-SIDH-master_M1 (or PQCrypto-SIDH-master_M2) for details.
------------------------------------------------------------------------------
Kaizhan Lin <linkzh5@mail2.sysu.edu.cn>