package scratch.kevin;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;

import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.gui.CPTListCellRenderer;
import org.opensha.commons.util.cpt.CPT;

public class CPTOptionViewer {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		ArrayList<CPT> cpts = GMT_CPT_Files.instances();
		
		JPanel panel = new JPanel();
		panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
		
		CPTListCellRenderer rend = new CPTListCellRenderer(null);
		for (CPT cpt : cpts) {
			panel.add(rend.buildComponent(cpt, 600, 30, Color.WHITE, Color.BLACK));
		}
		
		JFrame frame = new JFrame();
		frame.setContentPane(panel);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

}
