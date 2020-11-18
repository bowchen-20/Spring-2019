import java.io.IOException;

public class BlockDriver {
    public static void main(String[] args) throws IOException {
        Block block1 = new Block(20, 20, 400);
        block1.rdf("rdf.txt");
        block1.outConfig("config.txt");

        // block1.rdf("20000rdf.txt");
        // block1.outConfig("20000config.txt");

        block1.MD_simple(40000); 
        block1.outConfig("MD_equilibration.txt");
        block1.rdf("rdf_equilibration.txt");


        //---------T = 5K ----------------------
        // block1.setConfig("20000config.txt");
        // block1.MonteCarlo(50000, 5, "5K");
        // block1.outConfig("5K_config.txt");
        // block1.rdf("5K_rdf.txt");

        // block1.setConfig("5K_config.txt");
        // block1.MonteCarlo(50000, 10, "10K");
        // block1.outConfig("10K_config.txt");
        // block1.rdf("10K_rdf.txt");

        // block1.setConfig("10K_config.txt");
        // block1.MonteCarlo(50000, 15, "15K");
        // block1.outConfig("15K_config.txt");
        // block1.rdf("15K_rdf.txt");

        // block1.setConfig("15K_config.txt");
        // block1.MonteCarlo(50000, 20, "20K");
        // block1.outConfig("20K_config.txt");
        // block1.rdf("20K_rdf.txt");

        // block1.setConfig("20K_config.txt");
        // block1.MonteCarlo(50000, 25, "25K");
        // block1.outConfig("25K_config.txt");
        // block1.rdf("25K_rdf.txt");

        // block1.setConfig("25K_config.txt");
        // block1.MonteCarlo(50000, 30, "30K");
        // block1.outConfig("30K_config.txt");
        // block1.rdf("30K_rdf.txt");

        // block1.setConfig("30K_config.txt");
        // block1.MonteCarlo(50000, 35, "35K");
        // block1.outConfig("35K_config.txt");
        // block1.rdf("35K_rdf.txt");

        // block1.setConfig("35K_config.txt");
        // block1.MonteCarlo(50000, 40, "40K");
        // block1.outConfig("40K_config.txt");
        // block1.rdf("40K_rdf.txt");

        // block1.setConfig("40K_config.txt");
        // block1.MonteCarlo(50000, 45, "45K");
        // block1.outConfig("45K_config.txt");
        // block1.rdf("45K_rdf.txt");

        // block1.setConfig("45K_config.txt");
        // block1.MonteCarlo(50000, 50, "50K");
        // block1.outConfig("50K_config.txt");
        // block1.rdf("50K_rdf.txt");

        // block1.setConfig("50K_config.txt");
        // block1.MonteCarlo(50000, 55, "55K");
        // block1.outConfig("55K_config.txt");
        // block1.rdf("55K_rdf.txt");

        // block1.setConfig("55K_config.txt");
        // block1.MonteCarlo(50000, 60, "60K");
        // block1.outConfig("60K_config.txt");
        // block1.rdf("60K_rdf.txt");



        // //--- Now collect the equalibrium info from each state----
        // //--- Corresponding to each temperature ------------------

        // //------------T = 5K --------------------
        // block1.setConfig("5K_config.txt");
        // block1.CollectEqdata(10000, 5, "5K_Eq");

        // block1.setConfig("10K_config.txt");
        // block1.CollectEqdata(10000, 10, "10K_Eq");

        // block1.setConfig("15K_config.txt");
        // block1.CollectEqdata(10000, 15, "15K_Eq");

        // block1.setConfig("20K_config.txt");
        // block1.CollectEqdata(10000, 20, "20K_Eq");

        // block1.setConfig("25K_config.txt");
        // block1.CollectEqdata(10000, 25, "25K_Eq");

        // block1.setConfig("30K_config.txt");
        // block1.CollectEqdata(10000, 30, "30K_Eq");

        // block1.setConfig("35K_config.txt");
        // block1.CollectEqdata(10000, 35, "35K_Eq");

        // block1.setConfig("40K_config.txt");
        // block1.CollectEqdata(10000, 40, "40K_Eq");

        // block1.setConfig("45K_config.txt");
        // block1.CollectEqdata(10000, 45, "45K_Eq");

        // block1.setConfig("50K_config.txt");
        // block1.CollectEqdata(10000, 50, "50K_Eq");

        // block1.setConfig("55K_config.txt");
        // block1.CollectEqdata(10000, 55, "55K_Eq");

        // block1.setConfig("60K_config.txt");
        // block1.CollectEqdata(10000, 60, "60K_Eq");




    }
}
