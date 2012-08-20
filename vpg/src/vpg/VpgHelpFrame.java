/*
 * This file is part of the Visualization of Protein Ligand Graphs (VPLG) software package.
 *
 * Copyright Tim Schäfer 2012. VPLG is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */
package vpg;

/**
 *
 * @author ts
 */
public class VpgHelpFrame extends javax.swing.JFrame {

    /**
     * Creates new form VpgHelpFrame
     */
    public VpgHelpFrame() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jPanelContent = new javax.swing.JPanel();
        jLabelHelp = new javax.swing.JLabel();
        jScrollPaneHelpInfo = new javax.swing.JScrollPane();
        jTextAreaHelpInfo = new javax.swing.JTextArea();

        setTitle("VPG -- Help");

        jLabelHelp.setText("Getting help for VPG");

        jTextAreaHelpInfo.setEditable(false);
        jTextAreaHelpInfo.setBackground(new java.awt.Color(200, 200, 200));
        jTextAreaHelpInfo.setColumns(20);
        jTextAreaHelpInfo.setRows(5);
        jTextAreaHelpInfo.setText("Basic information on the VPLG programs can be found in\nthe doc/ directory of your VPLG installation.\n\nThis program has tooltips to help you, so just hover your \nmouse over an option you need help on.\n\nThe full documentation for VPLG is available online.");
        jScrollPaneHelpInfo.setViewportView(jTextAreaHelpInfo);

        javax.swing.GroupLayout jPanelContentLayout = new javax.swing.GroupLayout(jPanelContent);
        jPanelContent.setLayout(jPanelContentLayout);
        jPanelContentLayout.setHorizontalGroup(
            jPanelContentLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanelContentLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanelContentLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPaneHelpInfo, javax.swing.GroupLayout.DEFAULT_SIZE, 376, Short.MAX_VALUE)
                    .addGroup(jPanelContentLayout.createSequentialGroup()
                        .addComponent(jLabelHelp)
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
        jPanelContentLayout.setVerticalGroup(
            jPanelContentLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanelContentLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabelHelp)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPaneHelpInfo, javax.swing.GroupLayout.DEFAULT_SIZE, 255, Short.MAX_VALUE)
                .addContainerGap())
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanelContent, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanelContent, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(VpgHelpFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(VpgHelpFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(VpgHelpFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(VpgHelpFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new VpgHelpFrame().setVisible(true);
            }
        });
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel jLabelHelp;
    private javax.swing.JPanel jPanelContent;
    private javax.swing.JScrollPane jScrollPaneHelpInfo;
    private javax.swing.JTextArea jTextAreaHelpInfo;
    // End of variables declaration//GEN-END:variables
}
