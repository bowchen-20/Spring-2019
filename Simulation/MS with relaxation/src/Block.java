import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static java.lang.Math.*;

//Credits to Zhigao for helping me with the RDF
//https://matplotlib.org/1.3.1/users/usetex.html for details using latex

public class Block {
    private int w; // Width
    private int h; // Height
    private int n; // number of atoms

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

    // an atom is represented by an array
    private double[][] allAtoms;

    //a list of a nearby atoms of every single atom
    private List<Integer>[] neighborList;

    // compute Lenard Jones potential
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

    // compute the derivative of Lenard Jones potential
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
        allAtoms = new double[400][7];

        //set width and height to x,y
        w = x;
        h = y;
        this.n = n;

        for (int i = 0; i < n; i++)
        {
            double xpos, ypos;
            xpos = (i % w) * rMin; // default spacing between atoms to be rMin
            ypos = (i / h) * rMin;
            double[] atom = new double[7];
            atom[0] = xpos;
            atom[1] = ypos;
            atom[2] = atom[3] = atom[4] = atom[5] = atom[6] = 0.0;
            allAtoms[i] = atom;
        }

        // set the initial neighbor list
        getNeighborList(neighCut);

        // calculate the initial atomic stresses
        getAtomicStress();

        // calculate the initial Force on each atom
        calcForce();
    }

    ///------------------------------------------------------------
    /// getter functions
    ///------------------------------------------------------------
    int getN() {
        return allAtoms.length;
    }

    // return the distance between two atoms ----r_{ik}
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

    // return the distance in x axis between two atoms 
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

    // return the distance in y axis between two atoms 
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

    //------------get_neighbor_list---------------------------------
    // take a single argument r, will update neighbor_list of current configuration
    // -------------------------------------------------------------
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
                // r = neighCut
                if (getDistance(i, j) <= r)
                {
                    neighborList[i].add(j); // i's new neighbor is j
                    neighborList[j].add(i); // j's new neighbor is i
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
    // -------------------------------------------------------------
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

        // count each distance
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

        // last calculate g(r) and update the glist
        for (int i = 0; i < gList.size(); i++)
        {
            gList.set(i, gList.get(i) * 400 / omega / (2.0 * PI * deltaR * rList.get(i)));
            st = st + rList.get(i) + "   " + gList.get(i) + "\n";
        }

        // output rdf into a txt file for graphing
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
    // output the configuration of atoms for the plots
    // -------------------------------------------------------------
    void outConfig(String fileName) {
        String st = "";
        try {
            // record xpos, ypos and stresses of each atom into a txt file
            // file for latter plotting
            for (int i = 0; i < n; i++) {
                st = st + allAtoms[i][0] + " " + allAtoms[i][1] + " " + allAtoms[i][2] + " ";
                st = st + allAtoms[i][3] + " " + allAtoms[i][4] + "\n";
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

    //------------total_energy()------------------------------------
    // compute the total energy of a configuration by sum over the 
    // the potential energy of all atoms
    // -------------------------------------------------------------
    double totalEnergy() {
        double E = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < neighborList[i].size(); j++)
            {
                E += LJ(getDistance(i, neighborList[i].get(j)));
            }
        }
        // Since E_tot = 1/2 sum_{i,j} LJ(r_{i,j})
        E = E / 2.0;
        return E;
    }

    //------------calc_force() & calc_single_force()----------------
    //compute the force acting on each in a configuration 
    //
    //calc_single_force()
    //compute the net force on a single atom and stored at allAtoms[5]
    // and [6]
    //
    //calMaxForce() calculate the maximum force a current allAtoms
    //by sorting the force component
    // -------------------------------------------------------------
    double calcForce() {
        double f = 0;
        double curFx = 0;
        double curFy = 0;
        double curF = 0;
        // calc Force on each atom and find the F_max
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
            // update the list
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
    void calcSingleForce(int idx) {
        double curFx = 0;
        double curFy = 0;
        // calc Force on each atom and find the F_max
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

    //------------calc_total_stress_pressure()----------------------
    //compute the total stress sig_xx, simga_yy and simga_xy and
    //the hydrostatic pressure P
    //--------------------------------------------------------------
    List<Double> calcTotalStressPressure() {
        double sxx, syy, sxy, p;
        // P: total hydrostatic pressure:
        // for each atom i: p_i = -1/2*sum(sxx_i +syy_i)
        // P is just a sum over p_i
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

    //------------SD(lambda)---------------------------------
    //perform a steepest decent minimization, and for each distinct 
    //simulation, perform 10000 iterations by default
    // -------------------------------------------------------------
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
                        //first update the neighborlist
                        getNeighborList(neighCut);
                        //update the stresses
                        getAtomicStress();
                        //calc hydrostatic pressure
                        List<Double> sp = calcTotalStressPressure();
                        //calc total energy
                        double e = totalEnergy();
                        //calc maximum force
                        double f = calMaxForce();
                        //record the configuration data
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
        }
        catch(IOException e) {
            System.out.println("An error occurred while writing to file");
            e.printStackTrace();
        }
    }

}
