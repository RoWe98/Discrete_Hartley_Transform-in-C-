#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int N = 1024;
const float PI = 3.1416;

class DHT{
public:
	DHT(void); //构造函数
	~DHT(void);//析构函数
	inline void swap(float &a,float &b); //交换数据
	void bitrp(float xreal[], float ximag[], int n); //位反置转换
	void FFT(float xreal[], float ximag[], int n);   //快速傅立叶转换
	void FFT_test();     //测试数据
};

DHT::DHT(void)
{
	return;
}

DHT::~DHT(void)
{
	return;
}

inline void DHT::swap(float &a, float &b)
{
	float t;
	t = a;
	a = b;
	b = t;
}

void DHT::bitrp(float xreal[], float ximag[], int n)
{
	// 位反转置换 Bit-reversal Permutation
	int i, j, a, b, p;

	for (i = 1, p = 0; i < n; i *= 2)
	{
		p++;
	}
	for (i = 0; i < n; i++)
	{
		a = i;
		b = 0;
		for (j = 0; j < p; j++)
		{
			b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
			a >>= 1;        // a = a / 2;
		}
		if (b > i)
		{
			swap(xreal[i], xreal[b]);
			swap(ximag[i], ximag[b]);
		}
	}
}

void DHT::FFT(float xreal[], float ximag[], int n)
{
	// 快速傅立叶变换，将复数 x 变换后仍保存在 x 中，xreal, ximag 分别是 x 的实部和虚部
	float wreal[N / 2], wimag[N / 2], treal, timag, ureal, uimag, arg;
	int m, k, j, t, index1, index2;

	bitrp(xreal, ximag, n);

	// 计算 1 的前 n / 2 个 n 次方根的共轭复数 W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
	arg = -2 * PI / n;
	treal = cos(arg);
	timag = sin(arg);
	wreal[0] = 1.0;
	wimag[0] = 0.0;
	for (j = 1; j < n / 2; j++)
	{
		wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
		wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
	}

	for (m = 2; m <= n; m *= 2)
	{
		for (k = 0; k < n; k += m)
		{
			for (j = 0; j < m / 2; j++)
			{
				index1 = k + j;
				index2 = index1 + m / 2;
				t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
				treal = wreal[t] * xreal[index2] - wimag[t] * ximag[index2];
				timag = wreal[t] * ximag[index2] + wimag[t] * xreal[index2];
				ureal = xreal[index1];
				uimag = ximag[index1];
				xreal[index1] = ureal + treal;
				ximag[index1] = uimag + timag;
				xreal[index2] = ureal - treal;
				ximag[index2] = uimag - timag;
			}
		}
	}
}



void DHT::FFT_test()
{
	float xreal[N] = {}, ximag[N] = {};

	int n = 8;
	int i = 0;
	printf("请输入数据 格式(实部 虚部) 且长度等于8 : \n");
	for (i = 0; i < 8; i++)
	{
		scanf("%f%f",xreal+i,ximag+i);
	}

	n = i;    // 要求 n 为 2 的整数幂
	while (i > 1)
	{
		if (i % 2)
		{
			printf("%d is not a power of 2! ", n);
		}
		i /= 2;
	}

	FFT(xreal, ximag, n);
	printf("=====================================\n");
	printf("FFT:    i	      实部       虚部 \n");
	for (i = 0; i < n; i++)
	{
		printf("     %4d       %8.4f    %8.4f ", i+1, xreal[i], ximag[i]);
		printf("\n");
	}
	printf("===================================== \n");

	
	printf("DHT:    i	    结果 \n");
	for (i = 0; i < n; i++)
	{
		printf("     %4d        %8.4f ", i+1, xreal[i]-ximag[i]);
		printf("\n");
	}
	printf("=====================================\n ");
}

int main()
{
	DHT dht;
	dht.FFT_test();
	return 0;
}
