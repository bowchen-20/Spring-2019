import java.io.*;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;
import java.util.List;
import static java.lang.Math.*;
//discussed with Zhigao Han and sought help from Yudin Ai

public class Block {
    private int w; //Width
    private int h; //Height
    private int n; //number of atoms

    double eps = 0.010323;
    double sig = 3.405;
    double rMin = 3.822;
    double rTail = 7.0;
    double rCut = 7.5;
    double A = -0.0068102128;
    double B = -0.0055640876;
    double rdfMin = 0.0;
    double rdfMax = 12.0;
    double deltaR = 0.01;
    double neighCut = 7.5;
    double PI = 3.1415926;
    double k_b = 8.617*(Math.pow(10,-5));
    double delta_t = Math.pow(10,-14);
    double mass = 39.95 * 1.66 * Math.pow(10,-27);

    //an atom is represented by an array
    private double[][] allAtoms;

    //a list of a nearby atoms of every single atom
    private List<Integer>[] neighborList;

    //compute Lenard Jones potential
    double LJ(double r) {
        double result;
        if (r >= rCut) {
            result = 0;
        }
        else if (r >= rTail) {
            result = A * pow(r - rCut, 3) + B * pow(r - rCut, 2);
        }
        else {
            result = 4.0 * eps * (pow(sig / r, 12) - pow(sig / r, 6));
        }

        return result;
    }

    //compute the derivative of Lenard Jones potential
    double LJF(double r) {
        double result;
        if (r >= rCut) {
            result = 0;
        } else if (r >= rTail) {
            result = 3.0 * A * pow(r - rCut, 2) + 2.0 * B * (r - rCut);
        } else {
            result = (24.0 * eps / r) * (pow(sig / r, 6) - 2.0 * pow(sig / r, 12));
        }
        return result;
    }


    ///------------constructor--------------------------------------
    public Block(int x, int y, int n) {
        neighborList = new ArrayList[400];
        allAtoms = new double[400][9];

        //set width and height to x,y
        w = x;
        h = y;
        this.n = n;

        for (int i = 0; i < n; i++)
        {
            double xpos, ypos;
            xpos = (i % w) * rMin; //default spacing between atoms to be rMin
            ypos = (i / h) * rMin;
            double[] atom = new double[9];
            atom[0] = xpos;
            atom[1] = ypos;
            atom[2] = atom[3] = atom[4] = atom[5] = atom[6] = atom[7] = atom[8] = 0.0;

            allAtoms[i] = atom;
        }

        //set the initial neighbor list
        getNeighborList(neighCut);

        //calculate the initial atomic stresses
        getAtomicStress();

        //calculate the initial Force on each atom
        calcForce();
    }

    ///------------------------------------------------------------
    ///getter functions
    ///------------------------------------------------------------
    int getN() {
        return allAtoms.length;
    }

    //return the distance between two atoms ----r_{ik}
    double getDistance(int idx1, int idx2) {
        List<Double> coor1 = new ArrayList<>(2);
        coor1.add(allAtoms[idx1][0]);
        coor1.add(allAtoms[idx1][1]);
        List<Double> coor2 = new ArrayList<>(2);
        coor2.add(allAtoms[idx2][0]);
        coor2.add(allAtoms[idx2][1]);
        double dis;
        dis = sqrt(pow(coor1.get(0) - coor2.get(0), 2) + pow(coor1.get(1) - coor2.get(1), 2));
        return dis;
    }

    //return the distance in x axis between two atoms 
    double getXDistance(int idx1, int idx2) {
        List<Double> coor1 = new ArrayList<>(2);
        coor1.add(allAtoms[idx1][0]);
        coor1.add(allAtoms[idx1][1]);
        List<Double> coor2 = new ArrayList<>(2);
        coor2.add(allAtoms[idx2][0]);
        coor2.add(allAtoms[idx2][1]);
        double dis;
        dis = coor1.get(0) - coor2.get(0);
        return dis;
    }

    //return the distance in y axis between two atoms ----r_{ik}^{beta}
    double getYDistance(int idx1, int idx2) {
        List<Double> coor1 = new ArrayList<>(2);
        coor1.add(allAtoms[idx1][0]);
        coor1.add(allAtoms[idx1][1]);
        List<Double> coor2 = new ArrayList<>(2);
        coor2.add(allAtoms[idx2][0]);
        coor2.add(allAtoms[idx2][1]);
        double dis;
        dis = coor1.get(1) - coor2.get(1);
        return dis;
    }

    //------------getNeighborList---------------------------------
    //take a single argument r, will update neighbor_list of current configuration
    //-------------------------------------------------------------
    void getNeighborList(double r) {
        //create a new list

        for (int i = 0; i < neighborList.length; i++)
        {
            neighborList[i] = new ArrayList();
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                //r = neighCut
                if (getDistance(i, j) <= r)
                {
                    neighborList[i].add(j); //i's new neighbor is j
                    neighborList[j].add(i); //j's new neighbor is i
                }
            }
        }
    }

    //------------get_atomic_stress---------------------------------
    void getAtomicStress() {
        double omega = pow(rMin * (w - 1) * (h - 1) / n, 2);

        for (int i = 0; i < n; i++)
        {
            //set previous stress to 0
            allAtoms[i][2] = allAtoms[i][3] = allAtoms[i][4] = 0.0;

            for (int j = 0; j < neighborList[i].size(); j++)
            {
                double r = getDistance(i, neighborList[i].get(j));
                double rx = getXDistance(i, neighborList[i].get(j));
                double ry = getYDistance(i, neighborList[i].get(j));
                //sig_xx
                allAtoms[i][2] += LJF(r) * rx * rx / (r * omega);

                //sig_yy
                allAtoms[i][3] += LJF(r) * ry * ry / (r * omega);

                //sig_xy
                allAtoms[i][4] += LJF(r) * rx * ry / (r * omega);
            }
        }
    }

    //------------rdf() --------------------------------------------
    //Calculate the radial distribution function and store it in a text file
    //-------------------------------------------------------------
    void rdf(String fileName) {
        String st = "";
        double omega = pow(rMin * (w - 1) * (h - 1) / n, 2);
        List<Double> disList = new ArrayList<>();
        List<Double> rList = new ArrayList<>();
        List<Double> gList = new ArrayList<>();
        getNeighborList(rdfMax);

        for (int i = 0; i < 400; i++)
        {
            for (int j = 0; j < neighborList[i].size(); j++)
            {
                disList.add(getDistance(i, neighborList[i].get(j)));
            }
        }

        //dislist is sorted in increment order
        Collections.sort(disList);

        double d = disList.get(0);

        //count each distance
        int c = 0;
        for (int i = 0; i < disList.size(); i++)
        {
            if (abs(disList.get(i) - d) <= deltaR)
            {
                c++;
            }
            else
            {
                //rlist is a list of possible distance of all neighbors
                //sorted in increment order
                rList.add(d);
                //glist is a list of frequency regarding r
                gList.add(1.0 * c / disList.size()); //divide dislist to normalize it
                c = 0;
                d = disList.get(i);
            }

            if (i == disList.size() - 1)
            {
                rList.add(d);
                gList.add(1.0 * c / disList.size()); 
                c = 0;
                d = disList.get(i);
            }
        }

        //last calculate g(r) and update the glist
        for (int i = 0; i < gList.size(); i++)
        {
            gList.set(i, gList.get(i) * 400 / omega / (2.0 * PI * deltaR * rList.get(i)));
            st = st + rList.get(i) + "   " + gList.get(i) + "\n";
        }

        //output rdf into a txt file for graphing
        try {
            FileWriter myFileWriter = new FileWriter(fileName);
            myFileWriter.write(st);
            myFileWriter.close();
        }
        catch(IOException e) {
            System.out.println("An error occurred while writing to file");
            e.printStackTrace();
        }

        //reset the neighbor list to r = neighCut
        getNeighborList(neighCut);
    }

    //------------out_config()--------------------------------------
    //output the configuration of atoms for the plots
    //-------------------------------------------------------------
    void outConfig(String fileName) {
        String st = "";
        try {
            //record xpos, ypos and stresses of each atom into a txt file
            //file for latter plotting
            for (int i = 0; i < n; i++) {
                st = st + allAtoms[i][0] + " " + allAtoms[i][1] + " " + allAtoms[i][2] + " ";
                st = st + allAtoms[i][3] + " " + allAtoms[i][4] + " " + allAtoms[i][5] + " ";
                st = st + allAtoms[i][6] + " " + allAtoms[i][7] + " " + allAtoms[i][8] + "\n";
            }
            FileWriter myFileWriter = new FileWriter(fileName);
            myFileWriter.write(st);
            myFileWriter.close();
        }
        catch(IOException e) {
            System.out.println("An error occurred while writing to file");
            e.printStackTrace();
        }
    }

    void setConfig(String fileName) {
        try {
            FileInputStream fis = new FileInputStream(fileName);
            Scanner sc = new Scanner(fis);
            int i = 0; 
            while (sc.hasNextLine()) {
                String line = sc.nextLine();
                String[] data = line.split("\\s+");
                allAtoms[i][0] = Double.parseDouble(data[0]);
                allAtoms[i][1] = Double.parseDouble(data[1]);
                allAtoms[i][2] = Double.parseDouble(data[2]);
                allAtoms[i][3] = Double.parseDouble(data[3]);
                allAtoms[i][4] = Double.parseDouble(data[4]);
                allAtoms[i][5] = Double.parseDouble(data[5]);
                allAtoms[i][6] = Double.parseDouble(data[6]);
                allAtoms[i][7] = Double.parseDouble(data[7]);
                allAtoms[i][8] = Double.parseDouble(data[8]);
                i++;
            }
            sc.close();
        }
        catch(IOException e) {
            System.out.println("An error occurred while writing to file");
            e.printStackTrace();
        }
    }


    //------------total_energy()------------------------------------
    //compute the total energy of a configuration by sum over the 
    //the potential energy of all atoms
    //-------------------------------------------------------------
    double totalEnergy() {
        double E = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < neighborList[i].size(); j++)
            {
                E += LJ(getDistance(i, neighborList[i].get(j)));
            }
        }
        //Since E_tot = 1/2 sum_{i,j} LJ(r_{i,j})
        E = E / 2.0;
        return E;
    }

    //------------calc_force() & calc_single_force()----------------
    //compute the force acting on each in a configuration 
    //-------------------------------------------------------------
    double calcForce() {
        double f = 0;
        double curFx = 0;
        double curFy = 0;
        double curF = 0;
        //calc Force on each atom and find the F_max
        for (int i = 0; i < n; i++)
        {
            curFx = 0;
            curFy = 0;
            for (int j = 0; j < neighborList[i].size(); j++)
            {
                double r_ij = getDistance(i, neighborList[i].get(j));
                double r_ij_x = getXDistance(i, neighborList[i].get(j));
                double r_ij_y = getYDistance(i, neighborList[i].get(j));
                curFx += -LJF(r_ij) * (r_ij_x / r_ij);
                curFy += -LJF(r_ij) * (r_ij_y / r_ij);
            }
            //update the list
            allAtoms[i][5] = curFx;
            allAtoms[i][6] = curFy;
            curF = sqrt(pow(curFx, 2) + pow(curFy, 2));
            if (f < curF)
            {
                f = curF;
            }
        }
        return f;
    }

    void updatePosition() {
        for(int idx = 0; idx < allAtoms.length;idx++) {
        //x_position
        allAtoms[idx][0] = allAtoms[idx][0] + allAtoms[idx][7]*delta_t + delta_t*delta_t/(2*mass)* allAtoms[idx][5];

        //y_position
        allAtoms[idx][1] = allAtoms[idx][1] + allAtoms[idx][8]*delta_t + delta_t*delta_t/(2*mass)* allAtoms[idx][6];
        }
    }

    void updateVelocity(double[] F_prevx, double[] F_prevy) {
        for(int idx = 0; idx< allAtoms.length; idx++){
        //v_x
        allAtoms[idx][7] = allAtoms[idx][7] + delta_t/(2*mass) *(allAtoms[idx][5]+ F_prevx[idx]);

        //v_y
        allAtoms[idx][8] = allAtoms[idx][8] + delta_t/(2*mass) *(allAtoms[idx][6]+ F_prevy[idx]);
        }
    }
    
    void MD_simple(int it) throws IOException {
        String st = "";
        String temp_stress = "";
        List<Double> totalStress = new ArrayList<>(4);
        double[] presArray = new double[10];
        double[] tempArray = new double[10];
        double[] F_prevx = new double[400];
        double[] F_prevy = new double[400];
        int counter = 0;
        for (int n = 0; n<10; n++){
            for(int j = 0; j< allAtoms.length; j++){
                allAtoms[j][7] = 0;
                allAtoms[j][8] = 0;
            }
            double KE= 0;
            double KEsum = 0;
            double Psum=0; 
        double PE = totalEnergy();
        calcForce();

        int i = 0;
        while(i < it){
            for(int j = 0; j< allAtoms.length; j++){
                F_prevx[j] = allAtoms[j][5];
                F_prevy[j] = allAtoms[j][6];
            }
            // first get the new positions
            updatePosition();
            // we update the neighbor_list
            getNeighborList(neighCut);
            // find the new force
            calcForce();
            // velocity
            updateVelocity(F_prevx, F_prevy);

            //calculate new PE and KE 
            PE = totalEnergy();
            KE = 0; //recalc the KE
            for (int j = 0; j< allAtoms.length; j++){
                KE += 0.5 * mass* (allAtoms[j][7] * allAtoms[j][7] + 
                        allAtoms[j][8]*allAtoms[j][8]);
            }
            st = st + KE + " " + PE + " " + counter + "\n";
            i++;
            counter++;
            if(i>20000){ // assuming the block will be equilibrated after 20000 iterations
                KEsum +=KE;
                totalStress = calcTotalStressPressure();
                Psum += totalStress.get(3);
            }
        }
        KEsum /= (it - 20000); // calculate KE after equilibration 
        Psum /= (it - 20000); // calculate KE after equilibration 
        double temp = KEsum/(400*k_b);
        getAtomicStress(); 
        totalStress = calcTotalStressPressure();
        presArray[n] = totalStress.get(3);
        tempArray[n] = temp;
        temp_stress = temp_stress + temp + " " + totalStress.get(0) + " " + totalStress.get(0) + " "
                + totalStress.get(1) + " " + " " + totalStress.get(2) + " " + totalStress.get(3) + " " + Psum + "\n";
        }
        File file = new File("temp_stress.txt"); 
        FileWriter myFileWriter = new FileWriter(file);
        myFileWriter.write(temp_stress);
        myFileWriter.close();

        File file1 = new File("MD_Energy.txt");
        FileWriter myFileWriter1 = new FileWriter(file1);
        myFileWriter1.write(st);
        myFileWriter1.close();
    }

    void MD(int it, double Tobj, double Tinitial, String fileName) {
        // this for MD with temperature dependence not included because the T=0 version 
        // took too long to run and thus not sure of its validity. The biggest difference 
        // between MD_simple and MD is the inclusion scaling factor K. 

    }
    
    void calcSingleForce(int idx) {
        double curFx = 0;
        double curFy = 0;
        //calc Force on each atom and find the F_max
        for (int j = 0; j < neighborList[idx].size(); j++)
        {
            double r_ij = getDistance(idx, neighborList[idx].get(j));
            double r_ij_x = getXDistance(idx, neighborList[idx].get(j));
            double r_ij_y = getYDistance(idx, neighborList[idx].get(j));
            curFx += -LJF(r_ij) * (r_ij_x / r_ij);
            curFy += -LJF(r_ij) * (r_ij_y / r_ij);
        }
        allAtoms[idx][5] = curFx;
        allAtoms[idx][6] = curFy;
    }
    double calMaxForce() {
        double maxF = 0;
        double curFx = 0;
        double curFy = 0;
        double curF = 0;
        for (int i = 0; i < n; ++i)
        {
            curFx = allAtoms[i][5];
            curFy = allAtoms[i][6];
            curF = sqrt(pow(curFx, 2) + pow(curFy, 2));
            if (maxF < curF)
            {
                maxF = curF;
            }
        }
        return maxF;
    }

    //------------calcTotalStressPressure()----------------------
    //compute the total stress sig_xx, simga_yy and simga_xy and
    //the hydrostatic pressure P
    //--------------------------------------------------------------
    List<Double> calcTotalStressPressure() {
        double sxx, syy, sxy, p;
        //P: total hydrostatic pressure:
        //for each atom i: p_i = -1/2*sum(sxx_i +syy_i)
        //P is just a sum over p_i
        sxx = syy = sxy = p = 0;
        for (int i = 0; i < n; ++i)
        {
            sxx += allAtoms[i][2];
            syy += allAtoms[i][3];
            sxy += allAtoms[i][4];
            p += -1.0 / 2.0 * (allAtoms[i][2] + allAtoms[i][3]);
        }
        List<Double> sumStress = new ArrayList<>(4);
        sumStress.add(sxx);
        sumStress.add(syy);
        sumStress.add(sxy);
        sumStress.add(p);
        return sumStress;
    }

    void Move(double alpha,int selectIdx,double mcXpos,double mcYpos) {
        Random rand = new Random();     
        int randIdx = rand.nextInt(400);
        selectIdx = randIdx;

        double xzeta =  2.0 * Math.random() - 1.0;
        double yzeta =  2.0 * Math.random() - 1.0;

        mcXpos = allAtoms[randIdx][0];
        mcYpos = allAtoms[randIdx][1];

        allAtoms[randIdx][0] += alpha * xzeta;
        allAtoms[randIdx][1] += alpha * yzeta;
    }

    void moveAgain(int selectIdx, double mcXpos, double mcYpos){
        allAtoms[selectIdx][0] = mcXpos;
        allAtoms[selectIdx][1] = mcYpos;
    }

    void MonteCarlo(int steps, double T, String fileName) {
        String st = "";
        int it = 0;
        double PE = 0;
        double PEave = 0;
        int counter = 0;
        double alpha = 1;
        double sumalpha = 0;
        double probacceptave = 0;
        try {
            while(it<steps){
                PE = totalEnergy();
                double newPE, probaccept;
                double mcXpos = 0;
                double mcYpos = 0; 
                double prob;
                int selectIdx = 0;
                //move the atom
                Move(alpha,selectIdx,mcXpos,mcYpos);

                //update neighborlist and the PE
                getNeighborList(neighCut);
                newPE = totalEnergy();

                //Metropolis method;
                if(newPE>PE){
                    probaccept = exp(-(newPE-PE)/(k_b*T));
                    prob = Math.random();
                
                    if(prob > probaccept){
                        //do not accept if prob > probaccept
                        moveAgain(selectIdx,mcXpos,mcYpos); 
                    }
                    else{
                        //else accept move 
                        PE = newPE;
                    }
                }

                else{
                    probaccept = 1;
                    //accept the move
                    PE = newPE;
                }
                
                sumalpha += alpha;
                PEave += PE;
                probacceptave += probaccept;

                if(counter > (steps/500)){
                    sumalpha /= (steps/500);
                    PEave /= (steps/500);
                    probacceptave /= (steps/500);

                    //according to general aspect we want probacceptave to be 
                    //between 0.2 and 0.4, and if alpha exceeds a certain range 
                    //the complexity will be too big to handle
                    if(probacceptave > 0.35 && alpha < 5){
                        //increase alpha
                        alpha*=1.05;
                    }
                    else if (probacceptave <0.25){
                    //decrease alpha
                        alpha/=1.05;
                    }
                    counter = 0;
                    st = st + PEave + " " + probacceptave + " " + sumalpha + "\n";
                }
                counter++;
                it++;
            }
            getAtomicStress();
            FileWriter myFileWriter = new FileWriter(fileName);
            myFileWriter.write(st);
            myFileWriter.close();
        }
        catch(IOException e) {
            System.out.println("An error occurred while writing to file");
            e.printStackTrace();
        }
    }

    void CollectEqdata(int steps,double T,String fileName) {
        String st = "";
        int it = 0;
        double KE = 400*k_b*T;
        double PE = 0;
        double alpha = 1;
        double PEave = 0;
        double Pave = 0;
        int counter = 0;
        double sumalpha = 0;
        double probacceptave = 0;
        try{ 
            while(it<10000){
                PE = totalEnergy();
                double newPE,probaccept;
                double mcXpos = 0;
                double mcYpos = 0;
                double prob;
                int selectIdx = 0;
                //move the atom
                Move(alpha,selectIdx,mcXpos,mcYpos);

                //update neighborlist and the PE
                getNeighborList(neighCut);
                newPE = totalEnergy();

                //check the acceptance prob;
                if(newPE>PE){
                    probaccept = exp(-(newPE-PE)/(k_b*T));
                    prob = Math.random();
                
                    if(prob > probaccept){
                        //do not accept if prob > probaccept
                        moveAgain(selectIdx,mcXpos,mcYpos);
                    }
                    else{
                        //otherwise accept move 
                        PE = newPE;
                    }
                }
                else{
                    probaccept = 1;
                    //accept the move
                    PE = newPE;
                }

                sumalpha += alpha;
                probacceptave += probaccept;

                if(counter > (steps/100)){
                    sumalpha /= (steps/100);
                    probacceptave /= (steps/100);

                    //according to general aspect we want probacceptave to be
                    //between 0.2 and 0.4, and if alpha exceeds a certain range
                    //the complexity will be too big to handle
                    if(probacceptave > 0.35 && alpha < 5){
                        //increase alpha
                        alpha*=1.05;
                    }
                    else if (probacceptave <0.25){
                        //decrease alpha
                        alpha/=1.05;
                    }
                }
                PEave += PE;
                getAtomicStress();
                List<Double> stress = calcTotalStressPressure();
                Pave += stress.get(3);
                counter++;
                it++;
            }
            PEave /= steps;
            Pave /= steps;
            st = st + PEave + " " + KE + " " + Pave + "\n";
            FileWriter myFileWriter = new FileWriter(fileName);
            myFileWriter.write(st);
            myFileWriter.close();
        }
        catch(IOException e) {
            System.out.println("An error occurred while writing to file");
            e.printStackTrace();
        }
    }

    void sD(double lambda, int it, String fileName) {
        String st = "";
        int j = 0;
        int k = 0;
        try {
            while (j < it) {
                for (int i = 0; i < n; i++) {
                    calcSingleForce(i);
                    // Move each atom following the direction
                    // of it's force
                    allAtoms[i][0] += lambda * allAtoms[i][5];
                    allAtoms[i][1] += lambda * allAtoms[i][6];
                }

                if (it >= 1000) {
                    if (k == it / 1000) {
                        // first update the neighborlist
                        getNeighborList(neighCut);
                        // update the stresses
                        getAtomicStress();
                        // calc hydrostatic pressure
                        List<Double> sp = calcTotalStressPressure();
                        // calc total energy
                        double e = totalEnergy();
                        // calc maximum force
                        double f = calMaxForce();
                        // record the configuration data
                        st = st + j + "   " + e + "   " + sp.get(0) + "  " + sp.get(1) + "  ";
                        st = st + sp.get(2) + "  " + sp.get(3) + "  " + f + "\n";
                        k = 0;
                    }
                } else {
                    getNeighborList(neighCut);
                    getAtomicStress();
                    List<Double> sp = calcTotalStressPressure();
                    double e = totalEnergy();
                    double f = calMaxForce();
                    st = st + j + "   " + e + "   " + sp.get(0) + "  " + sp.get(1) + "  ";
                    st = st + sp.get(2) + "  " + sp.get(3) + "  " + f + "\n";
                }

                j++;
                k++;
            }
            FileWriter myFileWriter = new FileWriter(fileName);
            myFileWriter.write(st);
            myFileWriter.close();
        } catch (IOException e) {
            System.out.println("An error occurred while writing to file");
            e.printStackTrace();
        }
    }
}
