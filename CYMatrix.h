// CYMatrix.h: interface for the CCYMatrix class.
//
//////////////////////////////////////////////////////////////////////

/********************************************/
/*                                          */
/* ��� ������ �����ϴ� Ŭ����              */
/*       Ver 2.1,                           */
/*                                          */
/********************************************/

#ifndef __CCYMATRIX_CLASS_HEADER__
#define __CCYMATRIX_CLASS_HEADER__

#include "stdafx.h"
#include <stdio.h>

typedef double * PDOUBLE;

class CCYRect {
public:
	int left, top;
	int right, bottom;
	/////////////////////////////////
	CCYRect(){left=top=right=bottom=0;}
	CCYRect(int l, int t, int r, int b){left=l;top=t;right=r;bottom=b;}
	CCYRect(CCYRect &rect){*this=rect;}
	~CCYRect(){}
	int Width(){return right-left+1;}
	int Height(){return bottom-top+1;}	
	int operator==(CCYRect rect) {
		return (left==rect.left && top==rect.top && Width()==rect.Width() && Height()==rect.Height());
	}
	int operator!=(CCYRect rect) {
		return !(*this==rect);
	}
};

class CCYMatrix  
{
public:		
	/*
	int Load(char *filename);
	int Save(char *filename);
	void Load(FILE *fp);
	void Save(FILE *fp);
	*/

	CCYMatrix();
	CCYMatrix(int nRow, int nCol); // ��� ���� ũ�⸦ �ְ� ����� �ʱ�ȭ
	CCYMatrix(CCYMatrix &); // �ٸ� ��� �ν��Ͻ��� �ʱ�ȭ

	void operator=(CCYMatrix &); // = ������ 
	virtual ~CCYMatrix(); //
	void Free();
	
	void Create(int nRow, int nCol, double d=0); // �� ���·� ������ �ν��Ͻ� �ʱ�ȭ(����)
	int Row(){return _nRow;} // ���� ����� �� ũ�⸦ ��´�.
	int Col(){return _nCol;} // ���� ����� �� ũ�⸦ ��´�.
	CCYMatrix GetRow(int ndx); // ������ ���� ��´�.(���� ���·� ����)
	CCYMatrix GetCol(int ndx); // ������ ���� ��´�.(���� ���·� ����)
	double MaxAbs(); // ��� ���� ��� �� ���� ū ���� ��´�(����ġ)
	double Max(); // ���� ū ���� ��´�
	double MinAbs(); // ��� ���� ��� �� ���� ���� ���� ��´�(����ġ)
	double Min(); // ���� ���� ���� ��´�
	int IsError();
	
	CCYMatrix operator*(CCYMatrix &Op2); // ����� �� (����)
	CCYMatrix operator/(double Op2);	
	CCYMatrix operator+(CCYMatrix &Op2); // ����� �� (����)
	CCYMatrix operator-(CCYMatrix &Op2); // ����� �� (����)
	CCYMatrix operator^(int i); // Power

	double Det(); // ��Ľ�
	CCYMatrix  UTri(); // ��ﰢ���
	CCYMatrix Inv(); // ����� 
	CCYMatrix t(); // Transpose
	
	static void Copy(CCYMatrix &to, CCYMatrix &from, CCYRect toRect, int stx=0, int sty=0);
	static double inner_product(CCYMatrix &m1, CCYMatrix &m2); // inner product

	void Print(char *title = 0);
	int IsDataNull();

	inline double * operator[](unsigned long i){
		return m_pData[i];
	}

private:
	int _nCol; // ��� �� �� 
	int _nRow; // ��� �� ��
	PDOUBLE *m_pData; // ��� ��ü
};

#endif 
