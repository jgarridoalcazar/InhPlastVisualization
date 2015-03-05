#ifndef UEGORAND_H
#define UEGORAND_H

////////////////////////////////////////////////////////////
// $Id: uegorand.h,v 2.5 1998/03/17 23:14:52 jelasity Exp $
// uegorand.h
// declarations, macros etc related to psuedo-random number
// generation
////////////////////////////////////////////////////////////
// modification history:
//
////////////////////////////////////////////////////////////


#define MASK	2147483647	// 2^31 - 1
#define PRIME	65539		// 2^16 + 3
#define SCALE	( 1.0 / MASK )
#define RSEED	*(INI.Seed())


// ---- double value between 0 and 1, excluding 1.
#define UegoRand() \
	(( RSEED = ( ( RSEED * PRIME) & MASK) ) * SCALE )

// ---- long integer value between low and high inclusive
#define UegoLongRand( low, high ) \
	( (long) ( (long)(low) + ((long)(high)-(low)+1) * UegoRand()) )

// ---- double value between low and high, exluding high
#define UegoDoubleRand( low, high ) \
	( (low) + ((double)(high)-(low)) * UegoRand() )

#endif
