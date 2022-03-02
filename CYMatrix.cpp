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
		std::cout << "�޸� �Ҵ� ����! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
		exit(-1);
	}
	
	for(int row=0; row<this->_nRow; row++) {
		this->m_pData[row] = new double [this->_nCol];
		if(!this->m_pData[row]) {
			std::cout << "�޸� �Ҵ� ����! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
			exit(-1);
		}
	}
	
}

// ���� ������...
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
		std::cout << "�޸� �Ҵ� ����! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
		exit(-1);
	}

	for(int row=0; row<this->_nRow; row++) {
		this->m_pData[row] = new double [this->_nCol];
		if(!this->m_pData[row]) {
			std::cout << "�޸� �Ҵ� ����! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
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
		std::cout << "�޸� �Ҵ� ����! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
		exit(-1);
	}

	for(int row=0; row<this->_nRow; row++) {
		this->m_pData[row] = new double [this->_nCol];
		if(!this->m_pData[row]) {
			std::cout << "�޸� �Ҵ� ����! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
			exit(-1);
		}
		memcpy(this->m_pData[row], Src.m_pData[row], Src._nCol*sizeof(double));
	}	
}

// d�� �밢���� ���� (�� ������ ��� 0���� �ʱ�ȭ)
//   �����ϸ� �밢������ ������ ������ �ʱ�ȭ
void CCYMatrix::Create(int nRow, int nCol, double d)
{		
	if(!IsDataNull()) Free();

	this->_nRow = nRow; // 1
	this->_nCol = nCol; // 2

	this->m_pData = new PDOUBLE[this->_nRow];
	if(!this->m_pData) {
		std::cout << "�޸� �Ҵ� ����! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
		exit(-1);
	}

	for(int row=0; row<this->_nRow; row++) {
		this->m_pData[row] = new double [this->_nCol];
		if(!this->m_pData[row]) {
			std::cout << "�޸� �Ҵ� ����! : CCYMatrix::CCYMatrix(int nRowCount, int nColCount)" << std::endl;
			exit(-1);
		}

		for(int col=0; col<this->_nCol; col++) {
			if(row==col) (*this)[row][col] = d; // �밢����
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
/*   ����� �� ������ ����                */
/*                                        */
/******************************************/
CCYMatrix CCYMatrix::operator*(CCYMatrix &Op2)
{	
	CCYMatrix tmp;

	// �� ����� ���ҷ��� �� ����� ���� �� ����� ���� ���ƾ� ��
	if(this->Col() != Op2.Row()) {
		std::cout << "�� ����� ���� �� ����� ���� �ٸ�" << std::endl;
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
/*   ����� �� ������ ����                */
/*                                        */
/******************************************/
CCYMatrix CCYMatrix::operator+(CCYMatrix &Op2)
{
	CCYMatrix tmp;
	
	// ���ϱ� ������ �� ����� ��� ũ�Ⱑ ��Ȯ�� ���ƾ� ��
	if(!(this->Row()==Op2.Row() && this->Col()==Op2.Col())) {
		std::cout << "����� ������ ��ġ���� ����" << std::endl;
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
/*   ����� �� ������ ����                */
/*                                        */
/******************************************/
CCYMatrix CCYMatrix::operator-(CCYMatrix &Op2)
{
	CCYMatrix tmp;
	
	// ���� ������ �� ����� ��� ũ�Ⱑ ��Ȯ�� ���ƾ� ��
	if(!(this->Row()==Op2.Row() && this->Col()==Op2.Col())) {
		std::cout << "����� ������ ��ġ���� ����" << std::endl;
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

// ����� i�� ���� ����� �����ش� (power...)
// ��, i�� 1����...
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

// ������� ���Ѵ� ///////////////////
// (���콺 ���� �ҰŹ��� �̿�)
CCYMatrix CCYMatrix::Inv()
{
	CCYMatrix ai, work;
	if(this->Row() != this->Col()) {
		std::cout << "��������� �ƴ�" << std::endl;
		exit(-1);
		//return ai;
	}
	
	ai.Create(this->Row(), this->Col(), 1);	 // �밢 ������ 1�� �ʱ�ȭ
	work.Create(this->Row(), this->Col()*2, 1); // �밢 ������ 1�� �ʱ�ȭ

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

	for(i=0; i<work.Row(); i++) { // ���콺 ���� �ҰŹ�
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

// Transpose ���� ////////////////////////////////////
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

// ���� ���� ��� //////////////////////////////////
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

// from ����� Ư�� �κ��� to����� Ư�� �κп� �����Ѵ�
// ����1 ==> ���� ��, �� ��ȣ�� 0���� ����
// ����2 ==> CCYRect���� (left, top, right, bottom)������ �ε����� ����
// ------------------------------------------------------------
// ��) Copy(to, from, CCYRect(0,0,2,2), 0,0) : from����� (0��,0��)-(2��,2��)�� to����� (0��,0��)-(2��,2��)�� ����
// ��) Copy(to, from, CCYRect(0,0,2,2)) : ���� ����
// ��) Copy(to, from, CCYRect(3,0,5,2), 0,0) : from����� (0��,0��)-(2��,2��)�� to����� (0��,3��)-(2��,5��)�� ����
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


// �� �ﰢ���
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

// ��Ľ�
double CCYMatrix::Det()
{
	double v = 1.0;
	CCYMatrix m;
	m = this->UTri();

	for(int i=0; i<m.Col() && i<m.Row(); i++)
		v *= m[i][i];
	return v;
}

// ���� �� ��Ҹ� ��´�(�� �� ��ü�� ����)
CCYMatrix CCYMatrix::GetRow(int ndx)
{
	CCYMatrix tmp;
	tmp.Create(1,this->Col());
	for(int col=0; col<this->Col(); col++)
		tmp[0][col] = (*this)[ndx][col];

	return tmp;
}

// ���� �� ��Ҹ� ��´�(�� �� ��ü�� ����)
CCYMatrix CCYMatrix::GetCol(int ndx)
{
	CCYMatrix tmp;
	tmp.Create(this->Row(), 1);
	for(int row=0; row<this->Row(); row++)
		tmp[row][0] = (*this)[row][ndx];
	
	return tmp;
}

// ����ġ�� ���� �� ���� ū ���� ���Ѵ�.
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

// ����ġ�� ���� �� ���� ���� ���� ���Ѵ�.
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

// ���� ū ���� ��´�
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

// ���� ���� ���� ��´�
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