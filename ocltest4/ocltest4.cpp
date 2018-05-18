// ocltest4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "bn.h"

int main()
{
	BN_t bn1;

	unsigned int num1 = 0x4567;
	printf("num1 %d\t", num1);
	num1 -= 0x4500;
	printf("num1 %d\t", num1);
	num1 -= 0x78;
	if (num1 < 0)
		printf("num1 below zero\n");
	else
		printf("num1 still above zero\n");
	
	unsigned long long l1 = 0x1234FFFFFF;
	printf("test %llx\n", l1);
	//printf("test %d\n", 0x7FFFFFFF);

	unsigned char hex1[] = { 0x34, 0x56,0x67,0x78,0x9a,0xab,0xbc};
	BN_fromBuffer(&bn1, hex1, sizeof(hex1));
	printf("1 bignum length %d contains \n\t%08x\n\t%08x\n", bn1.length, bn1.words[0], bn1.words[1]);
	BN_fromBuffer(&bn1, hex1, sizeof(hex1));
	printf("2 bignum length %d contains \n\t%08x\n\t%08x\n", bn1.length, bn1.words[0], bn1.words[1]);

	BN_clone(&bn1, (BN_t *)&BN_nc);
	printf("3 bignum length %d contains \n\t%08x\n\t%08x\n\t%08x\n\t%08x\n\t%08x  \n", bn1.length, bn1.words[0], bn1.words[1], bn1.words[2], bn1.words[3], bn1.words[4]);


	BN_iuaddn(&bn1, 0xFF000000);
	printf("4 bignum length %d contains \n\t%08x\n\t%08x\n\t%08x\n\t%08x\n\t%08x  \n", bn1.length, bn1.words[0], bn1.words[1], bn1.words[2], bn1.words[3], bn1.words[4]);

    return 0;
}

/*
function ECPoint (x, y) {
	this.x = x
	this.y = y
	this.inf = false
}

ECPoint.fromPublicKey = function (publicKey) {
var first = publicKey[0]
var x
var y

if (publicKey.length === 33 && (first === 0x02 || first === 0x03)) {
x = BN.fromBuffer(publicKey.slice(1, 33))

// overflow
if (x.ucmp(BN.p) >= 0) return null

// create from X
y = x.redSqr().redMul(x).redIAdd7().redSqrt()
if (y === null) return null
if ((first === 0x03) !== y.isOdd()) y = y.redNeg()

return new ECPoint(x, y)
}

if (publicKey.length === 65 && (first === 0x04 || first === 0x06 || first === 0x07)) {
x = BN.fromBuffer(publicKey.slice(1, 33))
y = BN.fromBuffer(publicKey.slice(33, 65))

// overflow
if (x.ucmp(BN.p) >= 0 || y.ucmp(BN.p) >= 0) return null

// is odd flag
if ((first === 0x06 || first === 0x07) && y.isOdd() !== (first === 0x07)) return null

// x*x*x + 7 = y*y
if (x.redSqr().redMul(x).redIAdd7().ucmp(y.redSqr()) !== 0) return null

return new ECPoint(x, y)
}

return null
}

ECPoint.prototype.toPublicKey = function (compressed) {
	var x = this.x
	var y = this.y
	var publicKey

	if (compressed) {
		publicKey = Buffer.alloc(33)
		publicKey[0] = y.isOdd() ? 0x03 : 0x02
		x.toBuffer().copy(publicKey, 1)
	} else {
		publicKey = Buffer.alloc(65)
		publicKey[0] = 0x04
		x.toBuffer().copy(publicKey, 1)
		y.toBuffer().copy(publicKey, 33)
	}

	return publicKey
}

ECPoint.fromECJPoint = function (p) {
	if (p.inf) return new ECPoint(null, null)
	var zinv = p.z.redInvm()
	var zinv2 = zinv.redSqr()
	var ax = p.x.redMul(zinv2)
	var ay = p.y.redMul(zinv2).redMul(zinv)
	return new ECPoint(ax, ay)
}

ECPoint.prototype.toECJPoint = function () {
	if (this.inf) return new ECJPoint(null, null, null)
	return new ECJPoint(this.x, this.y, ECJPoint.one)
}

*/

/*
var g = require('./ecpointg')

ECPointG () {
	this.x = BN.fromBuffer(Buffer.from('79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798', 'hex'))
	this.y = BN.fromBuffer(Buffer.from('483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8', 'hex'))
	this.inf = false
	this._precompute()
}
ECPointG.prototype._precompute = function () {
	var ecpoint = new ECPoint(this.x, this.y)
	var dstep = 4
	var points = new Array(1 + Math.ceil(257 / dstep))
	var acc = points[0] = ecpoint
	for (var i = 1; i < points.length; ++i) {
		for (var j = 0; j < dstep; j++) acc = acc.dbl()
		points[i] = acc
	}

	this.precomputed = {
		naf: ecpoint._getNAFPoints(7),
		doubles: {
			step: dstep,
			points: points,
			negpoints: points.map(function (p) { return p.neg() })
		}
	}
}

ECPointG.prototype.mul = function (num) {
	var step = this.precomputed.doubles.step
	var points = this.precomputed.doubles.points
	var negpoints = this.precomputed.doubles.negpoints

	var naf = num.getNAF(1)
	var I = ((1 << (step + 1)) - (step % 2 === 0 ? 2 : 1)) / 3

	var repr = []
	for (var j = 0; j < naf.length; j += step) {
		var nafW = 0
		for (var k = j + step - 1; k >= j; k--) nafW = (nafW << 1) + naf[k]
		repr.push(nafW)
	}

	var a = new ECJPoint(null, null, null)
	var b = new ECJPoint(null, null, null)
	for (var i = I; i > 0; i--) {
		for (var jj = 0; jj < repr.length; jj++) {
			if (repr[jj] === i) {
				b = b.mixedAdd(points[jj])
			} else if (repr[jj] === -i) {
				b = b.mixedAdd(negpoints[jj])
			}
		}
		a = a.add(b)
	}

	return ECPoint.fromECJPoint(a)
}


g.mul(d).toPublicKey(compressed)



*/

/* 
function publicKeyTweakAdd(publicKey, tweak_bn, compressed) {
	var point = ECPoint.fromPublicKey(publicKey)
	return g.mul(tweak_bn).add(point).toPublicKey(compressed)
}
*/