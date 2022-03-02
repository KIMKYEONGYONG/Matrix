// CYMatrix.h: interface for the CCYMatrix class.
//
//////////////////////////////////////////////////////////////////////

/********************************************/
/*                                          */
/* 행렬 연산을 지원하는 클레스              */
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
	CCYMatrix(int nRow, int nCol); // 행과 열의 크기를 주고 행렬을 초기화
	CCYMatrix(CCYMatrix &); // 다른 행렬 인스턴스로 초기화

	void operator=(CCYMatrix &); // = 연산자 
	virtual ~CCYMatrix(); //
	void Free();
	
	void Create(int nRow, int nCol, double d=0); // 빈 상태로 시작한 인스턴스 초기화(생성)
	int Row(){return _nRow;} // 현재 행렬의 행 크기를 얻는다.
	int Col(){return _nCol;} // 현재 행렬의 열 크기를 얻는다.
	CCYMatrix GetRow(int ndx); // 지정된 행을 얻는다.(벡터 형태로 얻음)
	CCYMatrix GetCol(int ndx); // 지정된 행을 얻는다.(벡터 형태로 얻음)
	double MaxAbs(); // 행렬 내의 요소 중 가장 큰 값을 얻는다(절대치)
	double Max(); // 가장 큰 수를 얻는다
	double MinAbs(); // 행렬 내의 요소 중 가장 작은 값을 얻는다(절대치)
	double Min(); // 가장 작은 수를 얻는다
	int IsError();
	
	CCYMatrix operator*(CCYMatrix &Op2); // 행렬의 곱 (이항)
	CCYMatrix operator/(double Op2);	
	CCYMatrix operator+(CCYMatrix &Op2); // 행렬의 합 (이항)
	CCYMatrix operator-(CCYMatrix &Op2); // 행렬의 차 (이항)
	CCYMatrix operator^(int i); // Power

	double Det(); // 행렬식
	CCYMatrix  UTri(); // 상삼각행렬
	CCYMatrix Inv(); // 역행렬 
	CCYMatrix t(); // Transpose
	
	static void Copy(CCYMatrix &to, CCYMatrix &from, CCYRect toRect, int stx=0, int sty=0);
	static double inner_product(CCYMatrix &m1, CCYMatrix &m2); // inner product

	void Print(char *title = 0);
	int IsDataNull();

	inline double * operator[](unsigned long i){
		return m_pData[i];
	}

private:
	int _nCol; // 행렬 열 수 
	int _nRow; // 행렬 행 수
	PDOUBLE *m_pData; // 행렬 본체
};

#endif 
