////////////////////////////////////////////////////////////
// $Id: msgintf.cc,v 2.5 1998/03/17 23:14:55 jelasity Exp $
// msgintf.cc
// function for passing messages from the system to the user
////////////////////////////////////////////////////////////
// modification history:
//
////////////////////////////////////////////////////////////

#include "uego/uego.h"

#include <cstring>

short MsgLevel;

// -----------------------------------------------------------------------

void message(const char *msg, short level) {

	//auto start = std::chrono::high_resolution_clock::now();
	/*char white[15] = "\n\v\b\r\f";
	int pos;
	char tempMsg[2000];

	strcpy(tempMsg, msg);

	// ------- deleting whitespace from end --
	if (level > MSG_NOTHING) {
		pos = strcspn(tempMsg, white);
		tempMsg[pos] = 0;
		while (pos > 0 && (tempMsg[--pos] == ' ' || tempMsg[pos] == '\t'))
			tempMsg[pos] = 0;
	};

	switch (level) {
	case MSG_INFORMATION:
		if (MsgLevel >= MSG_INFORMATION)
			fprintf(stderr, "optipharm: -- %s\n", tempMsg);
		break;
	case MSG_ERROR:
		if (MsgLevel >= MSG_ERROR)
			fprintf(stderr, "optipharm: !! %s\n", tempMsg);
		break;
	default:
		if (MsgLevel >= MSG_NOTHING)
			fprintf(stderr, tempMsg);
	};
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;

	std::cout << "Elapsed time printf: " << elapsed.count() << " s\n";
*/
	//start = std::chrono::high_resolution_clock::now();
	switch (level) {
		case MSG_INFORMATION:
			clog << "optipharm: -- "<<msg<<"\n";
			break;
		case MSG_ERROR:
			cerr <<  "optipharm: !! "<<msg<<"\n";
			break;
		default:
			cout << msg<<"\n";
			break;
		};
	/*
	finish = std::chrono::high_resolution_clock::now();
		elapsed = finish - start;

		std::cout << "Elapsed time cout: " << elapsed.count() << " s\n";
		*/
}

//------------------------------------------------------------------------
/* setMsgLevel function
 *  it will set the MsgLevel variable for the correct value.
 */

void setMsgLevel(char level) {

	switch (level) {
	case '0':
		MsgLevel = 0;
		break;

	case '1':
		MsgLevel = 1;
		break;

	case '2':
		MsgLevel = 2;
		break;
	default:
		MsgLevel = 2; // default is 2 : display all msg
	}
}
