package scratch.kevin;

import javax.swing.JOptionPane;

public class AppNoLongerSupportedGUI {

	public static void main(String[] args) {
		String title = "Application No Longer Supported";
		String message = "This server mode application has been retired!"
				+ "\n"
				+ "\nNew standalone versions of each application are available"
				+ "\non our website: http://opensha.org/apps"
				+ "\n"
				+ "\nWill now exit.";
		
		JOptionPane.showMessageDialog(null, message, title, JOptionPane.ERROR_MESSAGE);
		System.exit(0);
	}

}
