package scratch.kevin;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;

import org.jdesktop.jxlayer.JXLayer;
import org.jdesktop.jxlayer.plaf.ext.LockableUI;

public class LockableUITest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		ActionListener gotClickedAL = new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				JButton source = (JButton)e.getSource();
				System.out.println("Clicked: "+source.getText());
			}
		};
		
		JFrame myFrame = new JFrame();
		JPanel panel = new JPanel();
		JButton button1 = new JButton("here's my button");
		button1.addActionListener(gotClickedAL);
		panel.add(button1);
		JPanel subPanel = new JPanel();
		JButton subButton1 = new JButton("button on a sub panel!");
		subPanel.add(subButton1);
		subButton1.addActionListener(gotClickedAL);
		JButton subPanelDisabledButton = new JButton("This subpanel button is already disabled.");
		subPanelDisabledButton.addActionListener(gotClickedAL);
		subPanelDisabledButton.setEnabled(false);
		subPanel.add(subPanelDisabledButton);
		panel.add(subPanel);
		
		final LockableUI lockUI = new LockableUI();
		JXLayer<JComponent> layer = new JXLayer<JComponent>(panel, lockUI);
		
		JPanel bigPanel = new JPanel(new BorderLayout());
		bigPanel.add(layer, BorderLayout.NORTH);
		JButton lockButton = new JButton("lock it!");
		lockButton.addActionListener(gotClickedAL);
		lockButton.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				lockUI.setLocked(!lockUI.isLocked());
			}
		});
		bigPanel.add(lockButton, BorderLayout.SOUTH);
		
		myFrame.setContentPane(bigPanel);
		myFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		myFrame.pack();
		myFrame.setVisible(true);
	}

}
