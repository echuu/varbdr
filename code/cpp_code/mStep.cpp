
#include "VarParam.h"
#include <stdio.h>

void VarParam::mStep() {


	printf("hello from mStep\n");

	this->m_0 = MAT_TYPE::Zero(D, K);



} // end mStep() function



int main(void) {
	printf("inside mStep file\n");
	return 0;
}

