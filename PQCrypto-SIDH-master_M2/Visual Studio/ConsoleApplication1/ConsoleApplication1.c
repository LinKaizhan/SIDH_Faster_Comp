// ConsoleApplication1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <stdio.h>
#include "../../src/config.h"
#include "../../src/P751/P751_internal.h"
#include "../../src/internal.h"

int main()
{
	const uint64_t B_basis_zero[8 * NWORDS64_FIELD] = { 0xF1A8C9ED7B96C4AB, 0x299429DA5178486E, 0xEF4926F20CD5C2F4, 0x683B2E2858B4716A, 0xDDA2FBCC3CAC3EEB, 0xEC055F9F3A600460,
													 0xD5A5A17A58C3848B, 0x4652D836F42EAED5, 0x2F2E71ED78B3A3B3, 0xA771C057180ADD1D, 0xC780A5D2D835F512, 0x114EA3B55AC1,
													 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
													 0x2E1EB8ED8C1C8C94, 0x6CFE456B25DBE01, 0x1EB54C3E8010F57A, 0x4B222D95FC81619D, 0xF99EBD204D501496, 0xC18348F9B629361,
													 0xC29E9A16BEDE6F96, 0x3B39F30163DAD41D, 0x807D3D1ECF2AC04E, 0xE088443F222A4988, 0x61B49A7524F1EA12, 0x41BF31133104,
													 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
													 0xE57361284693B54, 0xD66BD625AE87B791, 0x10B6D90DF32A3D0B, 0x97C4D1D7A74B8E95, 0x225D0433C353C114, 0x2AAA060C59FFB9F,
													 0xE46F50AF134F41D, 0x9442C2E31FC91DA1, 0xD920267A5E3844C3, 0xDDF0F4AD44A77A2A, 0x4691EACCBF84E753, 0x5E97318C9C5A,
													 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
													 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
													 0x2E1EB8ED8C1C8C94, 0x6CFE456B25DBE01, 0x1EB54C3E8010F57A, 0x4B222D95FC81619D, 0xF99EBD204D501496, 0xC18348F9B629361,
													 0xC29E9A16BEDE6F96, 0x3B39F30163DAD41D, 0x807D3D1ECF2AC04E, 0xE088443F222A4988, 0x61B49A7524F1EA12, 0x41BF31133104 };

	felm_t r;



	printf("B_basis_zero:\n");
	for (int k = 0; k < 8; k++)
	{
		from_mont((const digit_t*)&B_basis_zero + (k * NWORDS_FIELD), r);
		for (int i = 0; i < NWORDS_FIELD; i++)
		{
			printf("%016llx", r[NWORDS_FIELD - 1 - i]);
		}

		printf("\n");
		for (int i = 0; i < NWORDS_FIELD; i++)
		{
			printf("0x%llx,", r[i]);
		}
		printf("\n");
	}
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
