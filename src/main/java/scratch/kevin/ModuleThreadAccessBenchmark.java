package scratch.kevin;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.commons.util.modules.OpenSHA_Module;

import com.google.common.base.Stopwatch;

public class ModuleThreadAccessBenchmark {

	public static void main(String[] args) {
//		int numThreads = 20;
		int numThreads = 1;
		int accessesEach = 10000000;
		
		ModuleContainer<OpenSHA_Module> container = new ModuleContainer<>();
		
		container.addModule(new OpenSHA_Module() {
			
			@Override
			public String getName() {
				return "Module 1";
			}
		});
		container.addModule(new OpenSHA_Module() {
			
			@Override
			public String getName() {
				return "Module 2";
			}
		});
		container.addModule(new OpenSHA_Module() {
			
			@Override
			public String getName() {
				return "Module 3";
			}
		});
		container.addModule(new OpenSHA_Module() {
			
			@Override
			public String getName() {
				return "Module 4";
			}
		});
		container.addModule(new OpenSHA_Module() {
			
			@Override
			public String getName() {
				return "Module 5";
			}
		});
		container.addModule(new OpenSHA_Module() {
			
			@Override
			public String getName() {
				return "Module 6";
			}
		});
		container.addModule(new OpenSHA_Module() {
			
			@Override
			public String getName() {
				return "Module 7";
			}
		});
		container.addModule(new OpenSHA_Module() {
			
			@Override
			public String getName() {
				return "Module 8";
			}
		});
		container.addModule(new OpenSHA_Module() {
			
			@Override
			public String getName() {
				return "Module 9";
			}
		});
		container.addModule(new OpenSHA_Module() {
			
			@Override
			public String getName() {
				return "Module 10";
			}
		});
		
		System.out.println("Added "+container.getModules().size()+" modules");
		
		List<Class<? extends OpenSHA_Module>> moduleClasses = new ArrayList<>();
		for (OpenSHA_Module module : container.getModules())
			moduleClasses.add(module.getClass());
		
		List<AccessThread> threads = new ArrayList<>();
		for (int i=0; i<numThreads; i++)
			threads.add(new AccessThread(container, moduleClasses, accessesEach));
		
		Stopwatch watch = Stopwatch.createStarted();
		for (AccessThread thread : threads)
			thread.start();
		
		for (AccessThread thread : threads) {
			try {
				thread.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		watch.stop();
		double secs = (double)watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		System.out.println("Took "+(float)secs+" seconds");
	}
	
	private static class AccessThread extends Thread {
		
		private ModuleContainer<OpenSHA_Module> container;
		private List<Class<? extends OpenSHA_Module>> moduleClasses;
		private int numAccesses;

		public AccessThread(ModuleContainer<OpenSHA_Module> container,
				List<Class<? extends OpenSHA_Module>> moduleClasses, int numAccesses) {
			this.container = container;
			this.moduleClasses = moduleClasses;
			this.numAccesses = numAccesses;
		}
		
		@Override
		public void run() {
			Random r = new Random();
			for (int i=0; i<numAccesses; i++)
				container.getModule(moduleClasses.get(r.nextInt(moduleClasses.size())));
		}
		
	}

}
