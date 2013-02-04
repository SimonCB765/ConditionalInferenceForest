/**
 * 
 */
package tree;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Simon Bull
 *
 */
public class ProcessDataForGrowing
{

	/**
	 * The file containing the data from which the tree was grown.
	 */
	public String dataFileGrownFrom = "";
	/**
	 * 
	 */
	public Map<String, List<Object>> covariableData = new HashMap<String, List<Object>>();
	/**
	 * The number of observations in the data file.
	 */
	public int numberObservations;
	/**
	 * 
	 */
	public Map<String, List<Double>> responseData = new HashMap<String, List<Double>>();
	/**
	 * A set of the covariables that are used in the tree.
	 */
	public Set<String> covariablesGrownFrom;
	/**
	 * A mapping of the covariable names to their types.
	 * The type will be 'n', 'i', 'c' or 'x'.
	 */
	public Map<String, String> variableTypeMapping = new HashMap<String, String>();
	/**
	 * A mapping from the categorical variable names to the number of categories the variable has.
	 */
	public Map<String, Integer> categoricalVariableLevels = new HashMap<String, Integer>();

	public ProcessDataForGrowing()
	{
	}

	public ProcessDataForGrowing(String location)
	{
		try (BufferedReader reader = Files.newBufferedReader(Paths.get(location), StandardCharsets.UTF_8))
		{
			String line = reader.readLine();
			line = line.replaceAll("\n", "");
			String procDataVariables[] = line.split("\t");

			this.dataFileGrownFrom = procDataVariables[0];
			this.numberObservations = Integer.parseInt(procDataVariables[2]);
			this.variableTypeMapping = new HashMap<String, String>();
			if (!procDataVariables[5].equals(""))
			{
				String varMaping[] = procDataVariables[5].split(",");
				for (String s : varMaping)
				{
					String sSplit[] = s.split(";");
					this.variableTypeMapping.put(sSplit[0], sSplit[1]);
				}
			}
			this.covariableData = new HashMap<String, List<Object>>();
			if (!procDataVariables[1].equals(""))
			{
				String covData[] = procDataVariables[1].split(",");
				for (String s : covData)
				{
					String sSplit[] = s.split(";");
					List<Object> tmpList = new ArrayList<Object>();
					if (this.variableTypeMapping.get(sSplit[0]).equals("c"))
					{
						for (int i = 1; i < sSplit.length; i++)
						{
							tmpList.add(Integer.parseInt(sSplit[i]));
						}
					}
					else
					{
						for (int i = 1; i < sSplit.length; i++)
						{
							tmpList.add(Double.parseDouble(sSplit[i]));
						}
					}
					this.covariableData.put(sSplit[0], tmpList);
				}
			}
			this.responseData = new HashMap<String, List<Double>>();
			if (!procDataVariables[3].equals(""))
			{
				String respData[] = procDataVariables[3].split(",");
				for (String s : respData)
				{
					String sSplit[] = s.split(";");
					List<Double> tmpList = new ArrayList<Double>();
					for (int i = 1; i < sSplit.length; i++)
					{
						tmpList.add(Double.parseDouble(sSplit[i]));
					}
					this.responseData.put(sSplit[0], tmpList);
				}
			}
			this.covariablesGrownFrom = new HashSet<String>();
			if (!procDataVariables[4].equals(""))
			{
				String covGrown[] = procDataVariables[4].split(",");
				for (String s : covGrown)
				{
					this.covariablesGrownFrom.add(s);
				}
			}
			this.categoricalVariableLevels = new HashMap<String, Integer>();
			if (!procDataVariables[6].equals(""))
			{
				String covData[] = procDataVariables[6].split(",");
				for (String s : covData)
				{
					String sSplit[] = s.split(";");
					this.categoricalVariableLevels.put(sSplit[0], Integer.parseInt(sSplit[1]));
				}
			}
		}
		catch (Exception e)
		{
			System.err.println("Error: " + e.getMessage());
			System.exit(0);
		}
	}

	/**
	 * Copy constructor.
	 * 
	 * @param toCopy
	 */
	public ProcessDataForGrowing(ProcessDataForGrowing toCopy)
	{
		this.dataFileGrownFrom = toCopy.dataFileGrownFrom;
		this.covariableData = new HashMap<String, List<Object>>();
		this.covariableData.putAll(toCopy.covariableData);
		this.numberObservations = toCopy.numberObservations;
		this.responseData = new HashMap<String, List<Double>>();
		this.responseData.putAll(toCopy.responseData);
		this.covariablesGrownFrom = new HashSet<String>();
		this.covariablesGrownFrom.addAll(toCopy.covariablesGrownFrom);
		this.variableTypeMapping = new HashMap<String, String>();
		this.variableTypeMapping.putAll(toCopy.variableTypeMapping);
		this.categoricalVariableLevels = new HashMap<String, Integer>();
		this.categoricalVariableLevels.putAll(toCopy.categoricalVariableLevels);
	}

	/**
	 * Processes the data for growing.
	 *
	 * @return A class containing the mapping of covariable names to lists of their observed values, classifications to lists of their occurrences and the number of observations.
	 */
	public ProcessDataForGrowing(String data, TreeGrowthControl ctrl)
	{

		this.dataFileGrownFrom = data;
		ArrayList<ArrayList<Object>> vectorisedData = new ArrayList<ArrayList<Object>>();
		int responseColumn = 0;  // The variable indicating which column in the data file is the classifications.
		numberObservations = 0;
		String[] variableNames = null;  // The names of all the variables in the order they are in the file.
		String[] variableTypes = null;  // The types of all the variables in the order they are in the file.
		String[] variableSupports = null;  // The support of all the variables in the order they are in the file.
		Set<String> possibleClasses = new HashSet<String>();
		Path dataPath = Paths.get(dataFileGrownFrom);

		try (BufferedReader reader = Files.newBufferedReader(dataPath, StandardCharsets.UTF_8))
		{
			String line = null;

			// Record the names of the variables.
			line = reader.readLine();
			line = line.replaceAll("\n", "");
			variableNames = line.split("\t");

			// If the number of variables to exclude means that there are less than two variables remaining,
			// then alert the user and exit.
			if ((variableNames.length - ctrl.variablesToIgnore.size()) < 2)
			{
				System.out.println("You have masked out too many variables. You must leave at least one covariable and one response variable.");
				System.exit(0);
			}

			// Record the type of the variables in the file.
			// Currently assumes everything that is not a class is a numeric variable.
			line = reader.readLine();
			line = line.toLowerCase();
			line = line.replaceAll("\n", "");
			variableTypes = line.split("\t");

			// Record the specified support (number of categories) of the categorical variables.
			line = reader.readLine();
			line = line.toLowerCase();
			line = line.replaceAll("\n", "");
			variableSupports = line.split("\t", line.length() + 1);  // Explicit length so that the adjacent tabs do not get conflated (e.g. \t\t\t will generate [, , ] not [])
			// Determine whether only categorical variables have had their support defined.
			List<Integer> badlyDefinedSupport = new ArrayList<Integer>();
			for (int i = 0; i < variableTypes.length; i++)
			{
				if (!variableTypes[i].equals("c") && !variableSupports[i].equals("") && !variableTypes[i].equals("x"))
				{
					// If the variable is not categorical, and the support is defined.
					badlyDefinedSupport.add(i);
				}
			}
			if (!badlyDefinedSupport.isEmpty())
			{
				// If there are non-categorical variables with the support defined alert the user.
				System.out.print("The following set of columns are not for categorical variables, but have the number of categories defined : ");
				System.out.println(badlyDefinedSupport);
				System.exit(0);
			}

			int responseCount = 0;  // The number of variables that the user has indicated to be response variables.
			boolean isUnknownType = false;  // Whether or not an unknown type was encountered.
			for (int i = 0; i < variableTypes.length; i++)
			{
				if (variableTypes[i].equals("r"))
				{
					responseCount++;
					responseColumn = i;
				}
				else if (!variableTypes[i].equals("n") && !variableTypes[i].equals("i") && !variableTypes[i].equals("c") && !variableTypes[i].equals("x"))
				{
					isUnknownType = true;
				}
			}
			// If there was more or less than 1 class variable declared, then alert the user and exit.
			if (responseCount != 1)
			{
				System.out.format("You must define exactly one input variable to be the response. You defined %d.", responseCount);
				System.exit(0);
			}
			// If the response variable was masked out, then alert the user and exit.
			if (ctrl.variablesToIgnore.contains(variableNames[responseColumn]))
			{
				System.out.println("You have masked out the response variable.");
				System.exit(0);
			}
			// If an unknown type was encountered, then alert the user and exit.
			if (isUnknownType)
			{
				System.out.println("The only variable types permitted are 'r', 'c', 'n', 'i' and 'x'");
				System.exit(0);
			}

			// Turn the data into distinct columns, one for each variable that is not being ignored.
			for (int i = 0; i < variableNames.length; i++)
			{
				vectorisedData.add(new ArrayList<Object>());
			}
			while ((line = reader.readLine()) != null)
			{
				if (line.trim().length() == 0)
				{
					// If the line is made up of all whitespace, then ignore the line.
					continue;
				}
				numberObservations += 1;
				line = line.replaceAll("\n", "");
				String[] splitLine = line.split("\t");
				for (int i = 0; i < splitLine.length; i++)
				{
					if (i == responseColumn)
					{
						if (ctrl.isClassificationUsed)
						{
							possibleClasses.add(splitLine[i]);  // Record all the possible values that the class variable can take.
							vectorisedData.get(i).add(splitLine[i]);
						}
						else
						{
							vectorisedData.get(i).add(Double.parseDouble(splitLine[i]));
						}
					}
					else if (!ctrl.variablesToIgnore.contains(variableNames[i]))
					{
						// Record only the variables that are not marked as to be ignored (by either input or from the data file having an x for the variable type).
						if (variableTypes[i].equals("n") || variableTypes[i].equals("i"))
						{
							// If it is an integer or real valued covariable, then record it as a double.
							vectorisedData.get(i).add(Double.parseDouble(splitLine[i]));
						}
						else if (variableTypes[i].equals("c"))
						{
							// If it is a categorical covariable, then record it as an integer.
							vectorisedData.get(i).add(Integer.parseInt(splitLine[i]));
						}
					}
				}
			}

			// If there is only one class provided, then alert the user and exit.
			if (ctrl.isClassificationUsed && possibleClasses.size() < 2)
			{
				System.out.println("You have selected classification and the data has less than two classes.");
				System.exit(0);
			}
		}
		catch (IOException e)
		{
			// Caught an error while reading the file. Indicate this and exit.
			System.out.println("There was an error while processing the data file.\nGood Bye.");
			System.exit(0);
		}

		Map<String, String> tempCategoricalVariableLevels = new HashMap<String, String>();
		// Setup the mapping from variable names to variable types.
		// Also set up the temporary categorical variable category number mapping.
		for (int i = 0; i < variableNames.length; i++)
		{
			if (!ctrl.variablesToIgnore.contains(variableNames[i]))
			{
				variableTypeMapping.put(variableNames[i], variableTypes[i]);
				tempCategoricalVariableLevels.put(variableNames[i], variableSupports[i]);
			}
		}

		// Setup the mapping from variable name to the values of the variable for each observation.
		// Each entry in the Map will be an array of value. The nth value corresponds to the nth observation.
		// If you take the nth observation from each entry in the map, you have all data observed in the nth observation.
		for (int i = 0; i < variableNames.length; i++)
		{
			if (i != responseColumn && !variableTypes[i].equals("x"))
			{
				if (!ctrl.variablesToIgnore.contains(variableNames[i]))
				{
					covariableData.put(variableNames[i], vectorisedData.get(i));
				}
			}
		}
		this.covariablesGrownFrom = covariableData.keySet();

		if (ctrl.isClassificationUsed)
		{
			// Setup the matrix of the classification for each observation.
			// There is one column per class, and each class column has one entry for each observation.
			// If the value in a cell in the column i is 1, that means that the observation was classified as class i.
			// Each observation can only be classed as one class, and therefore for each observation there will be one column with a 1 and the rest will have a 0.
			// Example:
			// 		class1	class2	class3
			// obs1	0		1		0
			// obs2 1		0		0
			// ...
			// obsn 0		0		1
			List<String> classList = Arrays.asList(possibleClasses.toArray(new String[0]));  // A list of the possible classification values.
			List<Object> classificationValues = vectorisedData.get(responseColumn);  // A list of the classification for every observation.
			ArrayList<ArrayList<Double>> vectorisedClassData = new ArrayList<ArrayList<Double>>();  // A 2D array that will be converted into the Map strucutre to hold the final classification data.
			for (int i = 0; i < classList.size(); i++)
			{
				vectorisedClassData.add(new ArrayList<Double>());
			}
			for (int i = 0; i < numberObservations; i++)
			{
				Object classification = classificationValues.get(i);  // Get the actual classification.
				int classValueIndex = classList.indexOf(classification);  // Get the index of this class in the list of all possible classes.
				for (int j = 0; j < classList.size(); j++)
				{
					// For all possible values the class variable can take:
					if (j == classValueIndex)
					{
						// If the current class option is the classification for the current observation,
						// then set the entry to a 1. 
						vectorisedClassData.get(j).add(1.0);
					}
					else
					{
						// If the current class option is NOT the classification for the current observation,
						// then set the entry to a 0.
						vectorisedClassData.get(j).add(0.0);
					}
				}
			}
			for (int i = 0; i < classList.size(); i++)
			{
				String classValue = classList.get(i);
				int classValueIndex = classList.indexOf(classValue);
				responseData.put(classValue, vectorisedClassData.get(classValueIndex));
			}
		}
		else
		{
			// Classification is not used.
			List<Double> alteredResponseData = new ArrayList<Double>();
			List<Object> objectifiedResponseData = vectorisedData.get(responseColumn);
			for (int i = 0; i < numberObservations; i++)
			{
				alteredResponseData.add((Double) objectifiedResponseData.get(i));
			}
			responseData.put(variableNames[responseColumn], alteredResponseData);
		}

		// Setup the mapping from categorical variable names to the number of categories that the variable has.
		boolean isUserDefinedTooFewCategories = false;
		boolean isInorrectCategoryValue = false;
		List<String> badlyDefinedUserCategoryColumns = new ArrayList<String>();
		for (String i : variableTypeMapping.keySet())
		{
			if (variableTypeMapping.get(i).equals("c") && !ctrl.variablesToIgnore.contains(i))
			{
				// If the variable is categorical record the categories.
				List<Object> categoryData = covariableData.get(i);
				Set<Object> categories = new HashSet<Object>();
				for (Object j : categoryData)
				{
					categories.add(j);
				}
				// Check if the categories are ordered from 1 to the number of categories.
				// E.g. if there are three categories then they must be 1, 2 and 3.
				for (Object j : categories)
				{
					int valueOfJ = (Integer) j;
					if (valueOfJ < 1 || valueOfJ > categories.size())
					{
						// If the category value is less than 1 or greater than the
						// number of unique categories, then the categories are not defined in
						// the dataset as they should be.
						System.out.format("The categories for variable %s are not defined in the correct manner. They should be defined over a continuous range of integers, such that the smallest category value is 1 and the largest is the number of categories.\n", i);
						isInorrectCategoryValue = true;
						break;
					}
				}
				// If the use has specified the number of categories for this variable, make sure that it is not less than the number found.
				if (tempCategoricalVariableLevels.get(i).equals(""))
				{
					categoricalVariableLevels.put(i, categories.size());
				}
				else
				{
					int userDefinedCategoryNumber = Integer.parseInt(tempCategoricalVariableLevels.get(i));
					if (userDefinedCategoryNumber < categories.size())
					{
						// The user specified less categories than were found.
						isUserDefinedTooFewCategories = true;
						badlyDefinedUserCategoryColumns.add(i);
					}
					else
					{
						// The user specified no fewer categories than were found.
						categoricalVariableLevels.put(i, userDefinedCategoryNumber);
					}
				}
			}
		}
		if (isUserDefinedTooFewCategories)
		{
			System.out.println("The following set of variables had their number of categories defined, and the number supplied was fewer than were found in the dataset : ");
			System.out.println(badlyDefinedUserCategoryColumns);
			System.exit(0);
		}
		if (isInorrectCategoryValue)
		{
			System.exit(0);
		}

	}
	
	/**
	 * Method to save the control object to a file.
	 * 
	 * @param location The location where the file should be saved.
	 */
	public void save(String location)
	{
		try
		{
			FileWriter outputFile = new FileWriter(location);
			BufferedWriter outputWriter = new BufferedWriter(outputFile);
			outputWriter.write(this.dataFileGrownFrom + "\t");

			List<String> covarKeySet = new ArrayList<String>(this.covariableData.keySet());
			if (!covarKeySet.isEmpty())
			{
				outputWriter.write(covarKeySet.get(0));
				if (this.variableTypeMapping.get(covarKeySet.get(0)).equals("c"))
				{
					for (Object o : this.covariableData.get(covarKeySet.get(0)))
					{
						outputWriter.write(";");
						outputWriter.write(Integer.toString((Integer) o));
					}
				}
				else
				{
					for (Object o : this.covariableData.get(covarKeySet.get(0)))
					{
						outputWriter.write(";");
						outputWriter.write(Double.toString((Double) o));
					}
				}
				for (int i = 1; i < covarKeySet.size(); i++)
				{
					String currentCovar = covarKeySet.get(i);
					outputWriter.write("," + currentCovar);
					if (this.variableTypeMapping.get(currentCovar).equals("c"))
					{
						for (Object o : this.covariableData.get(currentCovar))
						{
							outputWriter.write(";");
							outputWriter.write(Integer.toString((Integer) o));
						}
					}
					else
					{
						for (Object o : this.covariableData.get(currentCovar))
						{
							outputWriter.write(";");
							outputWriter.write(Double.toString((Double) o));
						}
					}
				}
			}
			outputWriter.write("\t");

			outputWriter.write(Integer.toString(this.numberObservations) + "\t");

			List<String> responseKeySet = new ArrayList<String>(this.responseData.keySet());
			if (!responseKeySet.isEmpty())
			{
				outputWriter.write(responseKeySet.get(0));
				for (Double d : this.responseData.get(responseKeySet.get(0)))
				{
					outputWriter.write(";" + Double.toString(d));
				}
				for (int i = 1; i < responseKeySet.size(); i++)
				{
					outputWriter.write("," + responseKeySet.get(i));
					for (Double d : this.responseData.get(responseKeySet.get(i)))
					{
						outputWriter.write(";" + Double.toString(d));
					}
				}
			}
			outputWriter.write("\t");

			List<String> tmpCovars = new ArrayList<String>(this.covariablesGrownFrom);
			if (!tmpCovars.isEmpty())
			{
				outputWriter.write(tmpCovars.get(0));
				for (int i = 1; i < tmpCovars.size(); i++)
				{
					outputWriter.write("," + tmpCovars.get(i));
				}
			}
			outputWriter.write("\t");

			List<String> varTypeKeySet = new ArrayList<String>(this.variableTypeMapping.keySet());
			if (!varTypeKeySet.isEmpty())
			{
				outputWriter.write(varTypeKeySet.get(0) + ";" + this.variableTypeMapping.get(varTypeKeySet.get(0)));
				for (int i = 1; i < varTypeKeySet.size(); i++)
				{
					outputWriter.write("," + varTypeKeySet.get(i) + ";" + this.variableTypeMapping.get(varTypeKeySet.get(i)));
				}
			}
			outputWriter.write("\t");

			List<String> catLevelKeySet = new ArrayList<String>(this.categoricalVariableLevels.keySet());
			if (!catLevelKeySet.isEmpty())
			{
				outputWriter.write(catLevelKeySet.get(0) + ";" + Integer.toString(this.categoricalVariableLevels.get(catLevelKeySet.get(0))));
				for (int i = 1; i < catLevelKeySet.size(); i++)
				{
					outputWriter.write("," + catLevelKeySet.get(i) + ";" + Integer.toString(this.categoricalVariableLevels.get(catLevelKeySet.get(i))));
				}
			}
			outputWriter.close();
		}
		catch (Exception e)
		{
			System.err.println("Error: " + e.getMessage());
			System.exit(0);
		}
	}

}
