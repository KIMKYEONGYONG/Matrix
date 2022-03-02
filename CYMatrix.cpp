// CYMatrix.cpp: implementation of the CCYMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "CYMatrix.h"
#include <iostream>
#include <string.h> // for memcpy
#include <stdio.h>
#include <math.h> // for fabs ...
#include <process.h> // for exit(..)

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCYMatrix::CCYMatrix()
{
	this->_nRow = 0;    // 1 
	this->_nCol = 0;    // 2
	this->m_pData = NULL;

}

CCYMatrix::CCYMatrix(int nRowCount, int nColCount)
{
	this->_nRow = (nRowCount<0)?0:nRowCount; // 1
	this->_nCol = (nColCount<0)?0:nColCount; // 2

	this->m_pData = new PDOUBLE[this->_nRow];
	
	if(!this->m_pData) {
		std::cout << "메모리 할당 실패! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
		exit(-1);
	}
	
	for(int row=0; row<this->_nRow; row++) {
		this->m_pData[row] = new double [this->_nCol];
		if(!this->m_pData[row]) {
			std::cout << "메모리 할당 실패! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
			exit(-1);
		}
	}
	
}

// 복사 생성자...
CCYMatrix::CCYMatrix(CCYMatrix &Src)
{	
	this->m_pData = NULL;
	this->_nRow   = 0;    // 1
	this->_nCol   = 0;    // 2

	if(Src.IsDataNull()) return;

	this->_nRow = Src._nRow;    // 1
	this->_nCol = Src._nCol;    // 2

	this->m_pData = new PDOUBLE[this->_nRow];
	if(!this->m_pData) {
		std::cout << "메모리 할당 실패! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
		exit(-1);
	}

	for(int row=0; row<this->_nRow; row++) {
		this->m_pData[row] = new double [this->_nCol];
		if(!this->m_pData[row]) {
			std::cout << "메모리 할당 실패! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
			exit(-1);
		}
		memcpy(this->m_pData[row], Src.m_pData[row], Src._nCol*sizeof(double));
	}	

}

void CCYMatrix::operator=(CCYMatrix &Src)
{
	if(Src.IsDataNull()) return;
	if(!this->IsDataNull()) Free();

	this->_nRow = Src._nRow;    // 1
	this->_nCol = Src._nCol;    // 2
	
	this->m_pData = new PDOUBLE[this->_nRow];
	if(!this->m_pData) {
		std::cout << "메모리 할당 실패! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
		exit(-1);
	}

	for(int row=0; row<this->_nRow; row++) {
		this->m_pData[row] = new double [this->_nCol];
		if(!this->m_pData[row]) {
			std::cout << "메모리 할당 실패! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
			exit(-1);
		}
		memcpy(this->m_pData[row], Src.m_pData[row], Src._nCol*sizeof(double));
	}	
}

// d는 대각성분 지정 (미 지정시 모두 0으로 초기화)
//   지정하면 대각성분을 지정된 값으로 초기화
void CCYMatrix::Create(int nRow, int nCol, double d)
{		
	if(!IsDataNull()) Free();

	this->_nRow = nRow; // 1
	this->_nCol = nCol; // 2

	this->m_pData = new PDOUBLE[this->_nRow];
	if(!this->m_pData) {
		std::cout << "메모리 할당 실패! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
		exit(-1);
	}

	for(int row=0; row<this->_nRow; row++) {
		this->m_pData[row] = new double [this->_nCol];
		if(!this->m_pData[row]) {
			std::cout << "메모리 할당 실패! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
			exit(-1);
		}

		for(int col=0; col<this->_nCol; col++) {
			if(row==col) (*this)[row][col] = d; // 대각성분
			else (*this)[row][col] = 0.0;
		}
	}

}

CCYMatrix::~CCYMatrix()
{
	if(this->m_pData) {
		for(int row=0; row<Row(); row++)
			delete [] this->m_pData[row];
		delete [] this->m_pData;		
	}
	
	this->_nRow = 0;    // 1 
	this->_nCol = 0;    // 2
	this->m_pData = NULL;
}

void CCYMatrix::Free()
{
	if(this->m_pData) {
		for(int row=0; row<Row(); row++)
			delete [] this->m_pData[row];
		delete [] this->m_pData;		
	}
	
	this->_nRow = 0;    // 1 
	this->_nCol = 0;    // 2
	this->m_pData = NULL;
}

/******************************************/
/*                                        */
/*   행렬의 곱 연산을 수행                */
/*                                        */
/******************************************/
CCYMatrix CCYMatrix::operator*(CCYMatrix &Op2)
{	
	CCYMatrix tmp;

	// 두 행렬을 곱할려면 앞 행렬의 열과 뒷 행렬의 행이 같아야 함
	if(this->Col() != Op2.Row()) {
		std::cout << "앞 행렬의 열과 뒷 행렬의 행이 다름" << std::endl;
		exit(-1);
		//return tmp;
	}
		
	tmp.Create(this->Row(),Op2.Col());
	
	int targetRow = this->Row();
	int targetCol = Op2.Col();

	for(int row=0;row<targetRow;row++) {
		for(int col=0;col<targetCol;col++) {
			double mul_sum = 0;
			for(int entry=0;entry<Op2.Row();entry++)
				mul_sum += ((*this)[row][entry]*Op2[entry][col]);
			tmp[row][col] = mul_sum;
		}
	}

	return tmp;	
}


/******************************************/
/*                                        */
/*   행렬의 합 연산을 수행                */
/*                                        */
/******************************************/
CCYMatrix CCYMatrix::operator+(CCYMatrix &Op2)
{
	CCYMatrix tmp;
	
	// 더하기 연산은 두 행렬의 요소 크기가 정확히 같아야 함
	if(!(this->Row()==Op2.Row() && this->Col()==Op2.Col())) {
		std::cout << "행렬의 차원이 일치하지 않음" << std::endl;
		exit(-1);
		// return tmp;
	}
		
	tmp.Create(this->Row(),this->Col());
	
	int targetRow = this->Row();
	int targetCol = this->Col();

	for(int row=0;row<targetRow;row++) {
		for(int col=0;col<targetCol;col++) {
			tmp[row][col] = (*this)[row][col] + Op2[row][col];
		}
	}

	return tmp;	
}

/******************************************/
/*                                        */
/*   행렬의 차 연산을 수행                */
/*                                        */
/******************************************/
CCYMatrix CCYMatrix::operator-(CCYMatrix &Op2)
{
	CCYMatrix tmp;
	
	// 빼기 연산은 두 행렬의 요소 크기가 정확히 같아야 함
	if(!(this->Row()==Op2.Row() && this->Col()==Op2.Col())) {
		std::cout << "행렬의 차원이 일치하지 않음" << std::endl;
		exit(-1);
		//return tmp;
	}
		
	tmp.Create(this->Row(),this->Col());
	
	int targetRow = this->Row();
	int targetCol = this->Col();

	for(int row=0;row<targetRow;row++) {
		for(int col=0;col<targetCol;col++) {
			tmp[row][col] = (*this)[row][col] - Op2[row][col];
		}
	}

	return tmp;	
}

CCYMatrix CCYMatrix::operator/(double Op2)
{
	CCYMatrix tmp;
	
	tmp.Create(this->Row(),this->Col());
	
	int targetRow = this->Row();
	int targetCol = this->Col();

	for(int row=0;row<targetRow;row++) {
		for(int col=0;col<targetCol;col++) {
			tmp[row][col] = (*this)[row][col] / Op2;
		}
	}

	return tmp;	
}

// 행렬을 i번 곱한 결과를 돌려준다 (power...)
// 단, i는 1부터...
CCYMatrix CCYMatrix::operator^(int i)
{	
	if(i<=1) return *this;

	CCYMatrix tmp;
	tmp = *this;

	for(int t=1; t<i; t++) {
		tmp = tmp * (*this);
	}
	return tmp;
}

// 역행렬을 구한다 ///////////////////
// (가우스 조던 소거법을 이용)
CCYMatrix CCYMatrix::Inv()
{
	CCYMatrix ai, work;
	if(this->Row() != this->Col()) {
		std::cout << "정방행렬이 아님" << std::endl;
		exit(-1);
		//return ai;
	}
	
	ai.Create(this->Row(), this->Col(), 1);	 // 대각 성분을 1로 초기화
	work.Create(this->Row(), this->Col()*2, 1); // 대각 성분을 1로 초기화

	CCYMatrix::Copy(work,*this,CCYRect(0,0,this->Col()-1,this->Row()-1));
	CCYMatrix::Copy(work,ai,CCYRect(ai.Col(),0,ai.Col()*2-1,ai.Row()-1));

	int i, j, k, ctmp;
	double dtmp, pivot;

	for(i=0; i<work.Row()-1; i++) {
		ctmp = i;
		for(j=i+1; j<work.Row(); j++) {
			if(fabs(work[ctmp][i])<fabs(work[j][i])) ctmp = j;			
		}
		for(j=0; j<work.Col(); j++) {
			dtmp = work[i][j];
			work[i][j] = work[ctmp][j];
			work[ctmp][j] = dtmp;
		}
	}

	for(i=0; i<work.Row(); i++) { // 가우스 조던 소거법
		pivot = work[i][i];
		for(j=0; j<work.Col(); j++) {
			work[i][j] = work[i][j] / pivot;
		}
		for(j=0; j<work.Row(); j++) {
			if(i==j) continue;
			dtmp = work[j][i];
			for(k=i; k<work.Col(); k++)
				work[j][k] -= dtmp*work[i][k];
		}
	}

	Copy(ai, work, CCYRect(0,0,ai.Col()-1,ai.Row()-1),ai.Col(),0);

	return ai;
}

// Transpose 수행 ////////////////////////////////////
CCYMatrix CCYMatrix::t()
{
	CCYMatrix tmp;

	if(this->IsDataNull()) return tmp;

	tmp.Create(this->Col(), this->Row());

	for(int row=0; row<this->Row(); row++) {
		for(int col=0; col<this->Col(); col++) {			
			tmp[col][row] = (*this)[row][col];
		}
	}

	return tmp;
}

// 벡터 내적 계산 //////////////////////////////////
double CCYMatrix::inner_product(CCYMatrix &m1, CCYMatrix &m2)
{	
	if(m1.Col() != 1 || m2.Col() != 1) return -1.0;
	if(m1.Row() != m2.Row()) return -1.0;

	double mulSum = 0;
	for(int row=0; row<m1.Row(); row++) {
		mulSum += (m1[row][0]*m2[row][0]);
	}

	return mulSum;
}

// from 행렬의 특정 부분을 to행렬의 특정 부분에 복사한다
// 주의1 ==> 지정 행, 열 번호는 0부터 시작
// 주의2 ==> CCYRect에는 (left, top, right, bottom)순으로 인덱스를 지정
// ------------------------------------------------------------
// 예) Copy(to, from, CCYRect(0,0,2,2), 0,0) : from행렬의 (0행,0열)-(2행,2열)을 to행렬의 (0행,0열)-(2행,2열)로 복사
// 예) Copy(to, from, CCYRect(0,0,2,2)) : 위와 동일
// 예) Copy(to, from, CCYRect(3,0,5,2), 0,0) : from행렬의 (0형,0열)-(2행,2열)을 to행렬의 (0행,3열)-(2행,5열)로 복사
void CCYMatrix::Copy(CCYMatrix &to, CCYMatrix &from, CCYRect toRect, int stx, int sty)
{
	for(int toRow=toRect.top; toRow<=toRect.bottom; toRow++) {
		for(int toCol=toRect.left; toCol<=toRect.right; toCol++) {
			to[toRow][toCol] = from[sty+toRow-toRect.top][stx+toCol-toRect.left];
		}
	}
}

int CCYMatrix::IsDataNull()
{
	if(!this->m_pData) return 1;
	else return 0;
}

/////////////////////////////

void CCYMatrix::Print(char *title)
{
	if(this->IsDataNull()) return;

	if(title) {
		std::cout << title << std::endl;
	}

	for(int row=0; row<Row(); row++) {
		for(int col=0; col<Col(); col++) {
			printf("%+10.5lf ", (*this)[row][col]);
		}
		printf("\n");
	}

	std::cout << std::endl;
}


// 상 삼각행렬
CCYMatrix CCYMatrix::UTri()
{
	CCYMatrix at;
	at = *this;

	int i, j, k;
	double tmp, pivot;

	for(i=0; i<at.Col(); i++) {
		pivot = at[i][i];
		for(j=0; j<at.Col(); j++)
			at[i][j] /= pivot;
		for(j=i+1; j<at.Row(); j++) {
			tmp = at[j][i];
			for(k=i; k<at.Col(); k++)
				at[j][k] -= tmp*at[i][k];
		}
		for(j=0; j<at.Col(); j++)
			at[i][j] *= pivot;
	}

	return at;
}

// 행렬식
double CCYMatrix::Det()
{
	double v = 1.0;
	CCYMatrix m;
	m = this->UTri();

	for(int i=0; i<m.Col() && i<m.Row(); i++)
		v *= m[i][i];
	return v;
}

// 지정 행 요소를 얻는다(한 행 전체를 얻음)
CCYMatrix CCYMatrix::GetRow(int ndx)
{
	CCYMatrix tmp;
	tmp.Create(1,this->Col());
	for(int col=0; col<this->Col(); col++)
		tmp[0][col] = (*this)[ndx][col];

	return tmp;
}

// 지정 열 요소를 얻는다(한 열 전체를 얻음)
CCYMatrix CCYMatrix::GetCol(int ndx)
{
	CCYMatrix tmp;
	tmp.Create(this->Row(), 1);
	for(int row=0; row<this->Row(); row++)
		tmp[row][0] = (*this)[row][ndx];
	
	return tmp;
}

// 절대치로 했을 때 가장 큰 수를 구한다.
double  CCYMatrix::MaxAbs()
{
	double maxabs = fabs((*this)[0][0]);
	double theV = (*this)[0][0];
	for(int row=0; row<Row(); row++) {
		for(int col=0; col<Col(); col++) {
			if(maxabs<fabs((*this)[row][col])) {
				maxabs = fabs((*this)[row][col]);
				theV = (*this)[row][col];
			}
		}
	}
	return theV;
}

// 절대치로 했을 때 가장 작은 수를 구한다.
double  CCYMatrix::MinAbs()
{
	double minabs = fabs((*this)[0][0]);
	double theV = (*this)[0][0];
	for(int row=0; row<Row(); row++) {
		for(int col=0; col<Col(); col++) {
			if(minabs>fabs((*this)[row][col])) {
				minabs = fabs((*this)[row][col]);
				theV = (*this)[row][col];
			}
		}
	}
	return theV;
}

// 가장 큰 수를 얻는다
double  CCYMatrix::Max()
{
	double max = (*this)[0][0];
	for(int row=0; row<Row(); row++) {
		for(int col=0; col<Col(); col++) {
			if(max<(*this)[row][col]) {
				max = (*this)[row][col];
			}
		}
	}
	return max;
}

// 가장 작은 수를 얻는다
double  CCYMatrix::Min()
{
	double min = (*this)[0][0];
	for(int row=0; row<Row(); row++) {
		for(int col=0; col<Col(); col++) {
			if(min>(*this)[row][col]) {
				min = (*this)[row][col];
			}
		}
	}
	return min;
}

//////////////////////////////////////////////////////
/*
void CCYMatrix::Save(FILE *fp)
{
	fwrite(&this->_nRow, sizeof(this->_nRow), 1, fp);
	fwrite(&this->_nCol, sizeof(this->_nCol), 1, fp);
	
	for(int row=0; row<this->_nRow; row++) {
		fwrite(this->m_pData[row], sizeof(double)*this->_nCol, 1, fp);
	}
}

void CCYMatrix::Load(FILE *fp)
{
	if(!this->IsDataNull()) this->Free();

	fread(&this->_nRow, sizeof(this->_nRow), 1, fp);
	fread(&this->_nCol, sizeof(this->_nCol), 1, fp);

	this->m_pData = new PDOUBLE[this->_nRow];
	
	for(int row=0; row<this->_nRow; row++) {
		this->m_pData[row] = new double [this->_nCol];
		fread(this->m_pData[row], sizeof(double)*this->_nCol, 1, fp);
	}
}


int CCYMatrix::Save(char *filename)
{
	FILE *fp;
	if((fp=fopen(filename,"wb"))==NULL) return 0;

	Save(fp);
	
	fclose(fp);
	return 1; // no error
}

int CCYMatrix::Load(char *filename)
{
	FILE *fp;
	if((fp=fopen(filename,"rb"))==NULL) return 0;

	Load(fp);
	
	fclose(fp);
	return 1; // no error
}

*/