package scratch.kevin;

import mpi.MPI;
import mpi.MPIException;

public class MPJMemDebug {

	public static void main(String[] args) {
		try {
			MPI.Init(args);
			
			int me = MPI.COMM_WORLD.Rank();
			int size = MPI.COMM_WORLD.Size();
			
			Runtime rt = Runtime.getRuntime();
			long totalMB = rt.totalMemory() / 1024 / 1024;
			long maxMB = rt.maxMemory() / 1024 / 1024;
			System.out.println("Process "+me+"/"+size+": curMem="+totalMB+" MB, maxMem="+maxMB+" MB");
			
			MPI.Finalize();
		} catch (MPIException e) {
			MPI.COMM_WORLD.Abort(1);
		}
	}

}
