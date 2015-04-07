#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "uego.h"

short   MsgLevel;

// -----------------------------------------------------------------------

void    message( char *msg, short level ) {

        char    white[15] = "\n\v\b\r\f";
        int     pos;
        char    tempMsg[1000];

        time_t current_time;
        char  c_time_string[1000];

        strcpy(tempMsg,msg);

	// ------- deleting whitespace from end --
        if( level > MSG_NOTHING )
        {
		pos = strcspn( tempMsg, white );
                tempMsg[pos] = 0 ;
                while( pos>0 && (tempMsg[--pos]==' ' || tempMsg[pos]=='\t' ) )
                        tempMsg[pos] = 0;

        current_time = time(NULL);
        struct tm * p = localtime(&current_time);
        strftime(c_time_string, 1000, "%d/%m/%y -- %X", p);
        };

	switch( level )
	{
		case MSG_INFORMATION:
			if( MsgLevel >= MSG_INFORMATION )
				fprintf( stderr, "%s -- uego: -- %s\n", c_time_string, tempMsg );
			break;
		case MSG_ERROR:
			if( MsgLevel >= MSG_ERROR )
				fprintf( stderr, "%s -- uego: !! %s\n", c_time_string, tempMsg );
			break;
		default:
			if( MsgLevel >= MSG_NOTHING ) fprintf( stderr, tempMsg );
	};
}

//------------------------------------------------------------------------
/* setMsgLevel function
 *  it will set the MsgLevel variable for the correct value.
*/

void setMsgLevel( char level ) {

        switch(level)
        {
                case '0' :
                        MsgLevel = 0 ;
                        break ;

                case '1' :
                        MsgLevel = 1 ;
                        break ;

                case '2' :
                        MsgLevel = 2 ;
                        break ;
                default:
                        MsgLevel = 2 ; // default is 2 : display all msg
        }
}
