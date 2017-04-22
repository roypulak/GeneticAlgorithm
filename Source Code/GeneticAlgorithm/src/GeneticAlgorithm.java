import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class GeneticAlgorithm {
	private static int defaultGeneLength = 1000;
	private static byte[] bestGene = new byte[1000];
	private static byte[][] population = new byte[1000][defaultGeneLength];
	private static byte[][] tournamentPop = new byte[1000][defaultGeneLength];//used for individual selection before crossover

	private static int fittestIndex;// fittest index for population
	private static int tournamentFittestIndex;//fittest index for tournament Population
    private static final double uniformRate = 0.5;
    private static final double mutationRate = 0.015;
    private static final int tournamentSize = 5;
	private static final boolean elitism = true;
	
	private static int popSize = 200;
	private static int totalGeneration = 300;
	private static int totalNode = 7;
	
	private static int[][] G = new int [1000][1000];

	public static void main(String[] args) {
		
		double maxFitness,temp;
		int nthgen = 1;
		//extracting network from input file 
		processInput();
		
		// Create an initial population
		initialializePopulation(popSize, true);
		
		// Set a candidate solution
		maxFitness = getFittest(false);
		setBestIndividual();
		
		// Evolve our population until we reach an optimum solution
		for(int gen = 1; gen < totalGeneration; gen++)
		{
            System.out.println("Generation: " + gen + " Fitness: " + maxFitness);
            evolvePopulation();
            temp = getFittest(false);
            if(temp > maxFitness)
            {
            	maxFitness = temp;
            	setBestIndividual();
            	nthgen = gen+1;
            }	
		}
			
        System.out.println("Solution found!");
        System.out.println("Generation: " +nthgen);
        System.out.println("Genes:");
        for(int i = 0; i < defaultGeneLength; i++)
        System.out.print(bestGene[i]);
        writeFile();
	}
	
	// take best individual gene sequence
	private static void setBestIndividual()
	{
		for(int i = 0; i < defaultGeneLength; i++)
		{
			bestGene[i] = population[fittestIndex][i];
		}
	}
	
	  // Evolve a population
    public static void evolvePopulation() {
    	byte[] indv1 = new byte[defaultGeneLength];//used for crossover
    	byte[] indv2 = new byte[defaultGeneLength];//used for crossover
    	byte[] newIndiv = new byte [defaultGeneLength];//keeps new individual after crossover
 
        // Keep our best individual
        if (elitism) {
            for(int i = 0; i < defaultGeneLength; i++)
            {
            	population[0][i] = population[fittestIndex][i];
            }
        }

        // Crossover population
        int elitismOffset;
        if (elitism) {
            elitismOffset = 1;
        } else {
            elitismOffset = 0;
        }
        // Loop over the population size and create new individuals with
        // crossover
        for (int i = elitismOffset; i < popSize; i++) {
            int index1 = tournamentSelection(popSize);
            for(int j = 0; j < defaultGeneLength; j++)
            {
            	indv1[j] = tournamentPop[index1][j];
            }
            int index2 = tournamentSelection(popSize);
            for(int j = 0; j < defaultGeneLength; j++)
            {
            	indv2[j] = tournamentPop[index2][j];
            }

            newIndiv = crossover(indv1, indv2);
            for(int j = 0; j < defaultGeneLength; j++)
            {
            	population[i][j] = newIndiv[j];
            }
        }

        // Mutate population
        for (int i = elitismOffset; i < popSize; i++) {
            mutate(i);
        }
    }
    
    // Crossover individuals
    private static byte[] crossover(byte[] in1, byte[] in2) {
        byte[] newSol = new byte[defaultGeneLength];
        // Loop through genes
        for (int i = 0; i < in1.length; i++) {
            // Crossover
            if (Math.random() <= uniformRate) {
                newSol[i] = in1[i];
            } else {
                newSol[i] = in2[i];
            }
        }
        return newSol;
    }

    // Mutate an individual
    private static void mutate(int index) {
        // Loop through genes
        for (int i = 0; i < defaultGeneLength; i++) {
            if (Math.random() <= mutationRate) {
                // Create random gene
                byte gene = (byte) Math.round(Math.random());
                population[index][i] = gene;
            }
        }
    }

    // Select individuals for crossover
    private static int tournamentSelection(int pop) {
        // Create a tournament population
        //created tournament population in the above (globally) 
        // For each place in the tournament get a random individual
        for (int t = 0; t < tournamentSize; t++) {
            int randomId = (int) (Math.random() * pop); //here pop means popSize
            
            for(int i = 0 ; i < defaultGeneLength; i++)
            {
            	tournamentPop[t][i] = population[randomId][i];
            }
        }
        // Get the fittest among tournamentPop
        getFittest(true);
        return tournamentFittestIndex;// will return index
    }

    private static double getFittest(boolean isTournament) {
        double fitVal;
    	double temp;
    	
    	if(isTournament)// execute if called from tournamentSelection method
    	{          
        	fitVal = getFitness(0,true);
        	tournamentFittestIndex = 0;
        	
            // Loop through individuals to find fittest
            for (int i = 1; i < tournamentSize; i++) {
            	temp = getFitness(i,true);
                if (fitVal <= temp) {
                    fitVal = temp;
                    tournamentFittestIndex = i;
                }
            }
    	}
    	else
    	{
        	fitVal = getFitness(0,false);
        	fittestIndex = 0;
            // Loop through individuals to find fitVal
            for (int i = 1; i < popSize; i++) {
            	temp = getFitness(i,false);
                if (fitVal <= temp) {
                    fitVal = temp;
                    fittestIndex = i;
                }
            }
            
    	}
    	return fitVal;
    }

	public static void initialializePopulation(int populationSize, boolean initialise) {
		if (initialise) {
			// Loop and create individuals
			for (int i = 0; i < populationSize; i++) {
				for (int j = 0; j < defaultGeneLength; j++) {
					byte gene = (byte) Math.round(Math.random());
					population[i][j] = gene;
				}
			}
		}
	}

	static double getFitness(int index, boolean isTournament) {
		double sum = 0;
		double L_2;
		double Q_S;
		double F_S;
		//System.out.println("L2:"+L_2);
		L_2 = compute2L();
		// Loop through our individuals genes and compare them to our candidates
		if(isTournament)//if the method called from tournamentSelection() 
		{
			for (int i = 0;  i < totalNode; i++) {
				for(int j = 0; j < totalNode; j++)
				{
					sum += (G[i][j] - ( k(i) * k(j) ) / L_2 ) * ( tournamentPop[index][i] * tournamentPop[index][j] + (1-tournamentPop[index][i]) * (1-tournamentPop[index][j]));
				}
			}
		}
		else
		{
			for (int i = 0;  i < totalNode; i++) {
				for(int j = 0; j < totalNode; j++)
				{
					sum += (G[i][j] - ( k(i) * k(j) ) / L_2 ) * ( population[index][i] * population[index][j] + (1-population[index][i]) * (1-population[index][j]));
				}
			}
		}
		
		//Q_S = sum / L_2;
		Q_S = sum / L_2;
		//F_S = (1 + Q_S) / (double)2;
         F_S = Math.pow((1+Q_S), 4);
		return F_S;
	}
	
	//function that returns degree of a node
	static int k(int node)
	{
		int sum = 0;
		
		for(int i = 0; i < totalNode; i++)
		{
			if(G[node][i] == 1)
			{
				sum += G[node][i]; 
			}
				
		}
		
		return sum;
	}
	private static int compute2L()
	{
		int sum = 0;
		for(int i = 0; i < totalNode; i++)
		{
			for(int j = 0; j < totalNode; j++)
			{
				if(G[i][j] == 1)
				sum += G[i][j];
			}
		}
		return sum;
	}

	private static void processInput()
	{
		BufferedReader br = null;
		String[] parts = null;
		String line = null;
        int edgeCounter = 0;
        boolean isEdge = false;
        
        double tempN;
        //taking network file as input
		try {
			br = new BufferedReader(new FileReader("E:/Algorithm Codes/GeneticAlgorithm/networks/zachary_unwh.net"));
			while ((line = br.readLine()) != null) {
				if(line.contains("Vertices"))
				{
					parts = line.trim().split(" ");
					totalNode = Integer.parseInt(parts[1].trim());
					defaultGeneLength = totalNode;
					continue;
				}
				
				if(line.contains("Edges"))
				{
					isEdge = true;
					continue;
				}
				
				if(isEdge)
				{
					edgeCounter++;
					parts = line.trim().split(" ");
					//System.out.println(parts[0] + "  " + parts[1] + "  " + parts[2]);
					G[Integer.parseInt(parts[0]) - 1][Integer.parseInt(parts[1]) - 1] = 1;
					G[Integer.parseInt(parts[1]) - 1][Integer.parseInt(parts[0]) - 1] = 1;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}		
	}
	
	private static void writeFile()
	{
		try {

			//String content = "This is the content to write into file";

			File file = new File("E:/Algorithm Codes/GeneticAlgorithm/output/zachary_unwh.clu");

			// if file doesnt exists, then create it
			if (!file.exists()) {
				file.createNewFile();
			}

			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
		
			bw.write("*Vertices "+Integer.toString(totalNode)+"\n");
			for(int i = 0; i < defaultGeneLength; i++)
			{
				bw.write("\t "+Integer.toString(bestGene[i])+"\n");
				
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
