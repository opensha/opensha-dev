package scratch.kevin.simulators;

import java.io.IOException;

import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class QuickLoadTest {
	
	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance();
		
		catalog.loader().skipSlipsAndTimes().load();
	}

}
