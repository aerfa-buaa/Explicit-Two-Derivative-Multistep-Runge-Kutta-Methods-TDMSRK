// #pragma once
#include "xBurgers.hpp"

using namespace std;

void CSolver::Driver(void) {}

void CADSolver::Driver(void)
{

	InputControl_xy();

	MesherAndAllocation_xy();

	// Mesher_out();

	Initial_date();

	cout << " Solver Cell Number X Orientation = " << nCell << endl;
	cout << " Solver Cell Number Y Orientation = " << nCellYo << endl;
	omp_set_num_threads(NumberThreads);

	Boundary_uu();
	uu_to_cc();

	// out_date();
	if (Viscous_flag == 1)
		Viscous_mu();
	auto start = std::chrono::high_resolution_clock::now();
	auto start1 = std::chrono::high_resolution_clock::now();
	ExtIter = 0;
	TotalIter = 0;
	iSave = 0;
	TimeNow = 0.0;
	T_flag = 0.;
	while (ExtIter < StepMax)
	{

		ComputeDT();
		switch (TimeScheme)
		{

			// _________________________________
		case 33:
			RK33();
			break;
		case 54:
			RK54();
			break;
		case 65:
			RK65();
			break;
		case 76:
			RK76();
			break;
			// _________________________________

		case 223:
			TDMSRK223();
			break;
		case 224:
			TDMSRK224();
			break;

		case 235:
			TDMSRK235();
			break;

		case 326:
			TDMSRK326();
			break;
			// _________________________________
		case 231:
			TD23_New();
			break;
		case 241:
			TD24();
			break;

		case 351:
			TD35_New();
			break;

		case 461:
			TD46_New();
			break;

			// _________________________________

		default:
			cout << "TimeScheme is not avaliable for this Solver." << endl;
			break;
		}
		TimeNow += DT;

		// FinishTime = mydouble(clock()) / mydouble(CLOCKS_PER_SEC);
		// FinishTime = omp_get_wtime();
		// CurrentTotalTime = FinishTime - StartTime;
		if (ExtIter % show_n == 0)
		{
			cout << "Step = " << ExtIter << ";      Time = " << TimeNow << endl;
			auto stop1 = std::chrono::high_resolution_clock::now();
			auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1);
			std::cout << "Elapsed time: " << duration1.count() << " milliseconds" << std::endl;
		}
		ExtIter = ExtIter + 1;
		if (T_flag > 2.0)
			break;

		if (ska > 5)
		{
			cout << "Step = " << ExtIter << ";      Time = " << TimeNow << endl;
			cout << "NAN__NAN__NAN__NAN__NAN__NAN__= " << endl;
			break;
		}
	}
	cout << "Step = " << ExtIter - 1 << ";      Time = " << TimeNow << endl;
	auto stop_n = std::chrono::high_resolution_clock::now();
	// Calculate the elapsed time
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop_n - start);
	// Output the elapsed time in milliseconds
	std::cout << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;
	cout << "TimeScheme = " << TimeScheme << ";  CFL=" << CFL << endl;
	// #pragma omp barrier

	comput_err1();
	out_date();
	Delete_date();
}

void CADSolver::RK33(void)
{

	for (iRK_Step = 0; iRK_Step < 3; iRK_Step++)
	{

		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();
		TimeIntegration(iRK_Step);
		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::RK54(void)
{

	for (iRK_Step = 0; iRK_Step < 5; iRK_Step++)
	{

		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();

		TimeIntegration_RK54(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::RK65(void)
{

	for (iRK_Step = 0; iRK_Step < 6; iRK_Step++)
	{

		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();

		TimeIntegration_RK65(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}
void CADSolver::RK76(void)
{

	for (iRK_Step = 0; iRK_Step < 7; iRK_Step++)
	{

		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();

		TimeIntegration_RK76(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

// ________________________________

void CADSolver::TDMSRK223(void)
{

	if (ExtIter < 2)
	{
		// TD23();
		DT = DT * 0.5;
		TD23_New();
		if (ExtIter == 0)
			store_du_duu();
	}
	else
	{
		for (iRK_Step = 0; iRK_Step < 2; iRK_Step++)
		{

			Non_viscous_Splitting();

			// comput_du_LW_New();
			comput_du_LW_Old();
			scheme_du_dx_New();
			TimeIntegration_TDMSRK223_ALW(iRK_Step);

			Boundary_uu();
			uu_to_cc();
		}
	}
}

void CADSolver::TDMSRK224(void)
{

	if (ExtIter < 4)
	{
		DT = DT * 0.5;
		TD24();
		if (ExtIter == 0 || ExtIter == 2)
			store_du_duu();
	}
	else
	{
		for (iRK_Step = 0; iRK_Step < 2; iRK_Step++)
		{

			Non_viscous_Splitting();

			// comput_du_LW_New();
			comput_du_LW_Old();
			scheme_du_dx_New();
			TimeIntegration_TDMSRK224_ALW(iRK_Step);

			Boundary_uu();
			uu_to_cc();
		}
	}
}

void CADSolver::TDMSRK235(void)
{

	DT_w1 = DT_w2 = 1.0;
	if (ExtIter < 4)
	{
		DT = DT * 0.5;
		// TD23();
		TD35_New();
		if (ExtIter == 0 || ExtIter == 2 || ExtIter == 4)
			store_du_duu();
	}
	else
	{
		for (iRK_Step = 0; iRK_Step < 2; iRK_Step++)
		{

			Non_viscous_Splitting();

			// comput_du_LW_New();
			comput_du_LW_Old();
			scheme_du_dx_New();
			TimeIntegration_TDMSRK235_ALW(iRK_Step);
			// TimeIntegration_TDMSRK224_ALW(iRK_Step);
			// TimeIntegration_TDMS34_ALW(iRK_Step);

			Boundary_uu();
			uu_to_cc();
		}
	}
}

void CADSolver::TDMSRK326(void)
{
	if (ExtIter < 4)
	{
		DT = DT * 0.5;
		TD46_New();
		if (ExtIter == 0 || ExtIter == 2)
			store_du_duu();
	}
	else
	{
		for (iRK_Step = 0; iRK_Step < 3; iRK_Step++)
		{

			Non_viscous_Splitting();

			// comput_du_LW_New();
			comput_du_LW_Old();
			scheme_du_dx_New();
			TimeIntegration_TDMSRK326_ALW(iRK_Step);

			Boundary_uu();
			uu_to_cc();
		}
	}
}
// ________________________________

void CADSolver::TD23_New(void)
{

	for (iRK_Step = 0; iRK_Step < 2; iRK_Step++)
	{

		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();

		// comput_du_LW_New();
		comput_du_LW_Old();
		scheme_du_dx_New();
		if (Viscous_flag == 1)
			Viscous_Splitting_LW();

		TimeIntegration_TD23_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::TD24(void)
{

	for (iRK_Step = 0; iRK_Step < 2; iRK_Step++)
	{

		Non_viscous_Splitting();

		// comput_du_LW_New();
		comput_du_LW_Old();
		scheme_du_dx_New();
		TimeIntegration_TD24_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::TD35_New(void)
{
	for (iRK_Step = 0; iRK_Step < 3; iRK_Step++)
	{
		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();

		// comput_du_LW_New();
		comput_du_LW_Old();
		scheme_du_dx_New();
		if (Viscous_flag == 1)
			Viscous_Splitting_LW();

		TimeIntegration_TD35_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}

void CADSolver::TD46_New(void)
{
	for (iRK_Step = 0; iRK_Step < 4; iRK_Step++)
	{
		Non_viscous_Splitting();
		if (Viscous_flag == 1)
			Viscous_Splitting();

		// comput_du_LW_New();
		comput_du_LW_Old();
		scheme_du_dx_New();
		if (Viscous_flag == 1)
			Viscous_Splitting_LW();

		TimeIntegration_TD46_ALW(iRK_Step);

		Boundary_uu();
		uu_to_cc();
	}
}
//________________________________
