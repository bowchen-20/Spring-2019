public class BlockDriver {
    public static void main(String[] args) {
        Block block1 = new Block(20, 20, 400);
        // output the original configuration and rdf
        block1.rdf("rdf.txt");
        block1.outConfig("config.txt");

        // run lambda at 0.5, 1, 1.5, 5, 10 respectively

        // ---------------lambda = 0.5 -------------
        block1.sD(0.5, 10000, "l=0p5_vsN.txt");
        block1.outConfig("l=0p5_config.txt");
        block1.rdf("l=0p5_rdf.txt");

        // ---------------lambda = 1 ---------------
        Block block2 = new Block(20, 20, 400);
        block2.sD(1, 10000, "l=1_vsN.txt");
        block2.outConfig("l=1_config.txt");
        block2.rdf("l=1_rdf.txt");

        // ---------------lambda = 1.5 -------------
        Block block3 = new Block(20, 20, 400);
        block3.sD(1.5, 10000, "l=1p5_vsN.txt");
        block3.outConfig("l=1p5_config.txt");
        block3.rdf("l=1p5_rdf.txt");

        // ---------------lambda = 5 ---------------
        Block block4 = new Block(20, 20, 400);
        block4.sD(5, 10000, "l=5_vsN.txt");
        block4.outConfig("l=5_config.txt");
        block4.rdf("l=5_rdf.txt");

        // ---------------lambda = 10 ---------------
        Block block5 = new Block(20, 20, 400);
        block5.sD(10, 10000, "l=10_vsN.txt");
        block5.outConfig("l=10_config.txt");
        block5.rdf("l=10_rdf.txt");

        // --------study the evolution -------------
        // ---------------lambda = 1.5 -------------
        Block block6 = new Block(20, 20, 400);

        block6.sD(1.5, 100, "100vsN.txt");
        block6.outConfig("100config.txt");
        block6.rdf("100rdf.txt");

        block6.sD(1.5, 500, "500vsN.txt");
        block6.outConfig("500config.txt");
        block6.rdf("500rdf.txt");

        block6.sD(1.5, 1000, "1000vsN.txt");
        block6.outConfig("1000config.txt");
        block6.rdf("1000rdf.txt");

        block6.sD(1.5, 5000, "5000vsN.txt");
        block6.outConfig("5000config.txt");
        block6.rdf("5000rdf.txt");

        block6.sD(1.5, 10000, "10000vsN.txt");
        block6.outConfig("10000config.txt");
        block6.rdf("10000rdf.txt");

        block6.sD(1.5, 20000, "20000vsN.txt");
        block6.outConfig("20000config.txt");
        block6.rdf("20000rdf.txt");

        block6.sD(1.5, 50000, "50000vsN.txt");
        block6.outConfig("50000config.txt");
        block6.rdf("50000rdf.txt");
    }
}
