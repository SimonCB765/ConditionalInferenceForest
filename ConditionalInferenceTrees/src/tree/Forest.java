/**
 * 
 */
package tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * @author Simon Bull
 *
 */
public class Forest
{

	/**
	 * A list of the trees in the forest.
	 */
	public List<ConditionalInferenceTree> forest = new ArrayList<ConditionalInferenceTree>();
	/**
	 * The oob observations ordered such that the ith oob list is for the ith observation.
	 */
	public List<List<Integer>> oobObservations = new ArrayList<List<Integer>>();
	/**
	 * The oob error estimate.
	 */
	public double oobErrorEstimate = 0.0;
	/**
	 * The file containing the data that the forest was grown from.
	 */
	public String dataFileGrownFrom = "";
	/**
	 * The object recording the control parameters for the forest and its trees.
	 */
	public TreeGrowthControl ctrl;
	/**
	 * 
	 */
	public ProcessDataForGrowing processedData;

	public Forest(String dataForGrowing)
	{
		this.ctrl = new TreeGrowthControl();
		growForest(dataForGrowing, null);
	}

	public Forest(String dataForGrowing, TreeGrowthControl ctrl)
	{
		this.ctrl = ctrl;
		growForest(dataForGrowing, null);
	}

	public Forest(String dataForGrowing, Double[] weights)
	{
		this.ctrl = new TreeGrowthControl();
		growForest(dataForGrowing, weights);
	}

	public Forest(String dataForGrowing, TreeGrowthControl ctrl, Double[] weights)
	{
		this.ctrl = ctrl;
		growForest(dataForGrowing, weights);
	}

	void growForest(String dataForGrowing, Double[] potentialWeights)
	{

		ProcessDataForGrowing processedData = new ProcessDataForGrowing(dataForGrowing, this.ctrl);
		this.processedData = processedData;

		// Setup the record of oobObservations, one record for each observation.
		for (int i = 0; i < processedData.numberObservations; i++)
		{
			oobObservations.add(new ArrayList<Integer>());
		}

		// Generate the default weightings.
		Double[] weights;
		if (potentialWeights == null)
		{
			weights = new Double[processedData.numberObservations];
			Arrays.fill(weights, 1.0);
		}
		else
		{
			weights = potentialWeights;
		}

		// Setup the observation selection variables.
		List<Integer> observations = new ArrayList<Integer>();
		for (int i = 0; i < processedData.numberObservations; i++)
		{
			observations.add(i);
		}
		int numberObservationsToSelect = 0;
		if (!ctrl.isReplacementUsed)
		{
			numberObservationsToSelect = (int) Math.floor(ctrl.selectionFraction * processedData.numberObservations);
		}
		else
		{
			numberObservationsToSelect = processedData.numberObservations;
		}

		for (int i = 0; i < ctrl.numberOfTreesToGrow; i++)
		{
			// Randomly determine the observations used for growing this tree.
			Double[] weightsForThisTree = new Double[processedData.numberObservations];
			Arrays.fill(weightsForThisTree, 0.0);
			if (!ctrl.isReplacementUsed)
			{
				// Selection of observations when replacement is not used involves setting the weights of
				// a fraction of the variables to their original value, while the remainder get 0.
				Collections.shuffle(observations);
				for (int j = 0; j < numberObservationsToSelect; j++)
				{
					weightsForThisTree[observations.get(j)] = weights[observations.get(j)];
				}
			}
			else
			{
				// Selection of observations when replacement is used involves initialising the weights of
				// all variables to 0, and incrementing the weight of any variables chosen.
				Random randomObservation = new Random();
				int selectedObservation;
				for (int j = 0; j < numberObservationsToSelect; j++)
				{
					selectedObservation = randomObservation.nextInt(processedData.numberObservations);
					weightsForThisTree[selectedObservation] += weights[selectedObservation];
				}
			}

			// Update the list of which trees an observation is oob on.
			for (int j = 0; j < processedData.numberObservations; j++)
			{
				if (weightsForThisTree[j] == 0 && weights[j] != 0)
				{
					// If the jth observation has a weight of 0, then it is oob for this tree (the ith tree).
					// Unless the weight for the observation was originally set to 0. In this case the observation is
					// not to be used, and therefore not oob.
					oobObservations.get(j).add(i);
				}
			}

			// Grow this tree from the chosen observations.
			this.forest.add(new ConditionalInferenceTree(processedData, this.ctrl, weightsForThisTree));
		}

		// Calculate the oob error.
		double cumulativeErrors = 0.0;  // The number of errors observed over all the oob observations.
		Set<Integer> observationsThatAreOob = new HashSet<Integer>();
		for (int i = 0; i < processedData.numberObservations; i++)
		{
			List<Integer> observationsToPredict = new ArrayList<Integer>();
			observationsToPredict.add(i);
			if (!oobObservations.get(i).isEmpty())
			{
				// If this observation is oob at least once, then predict it on all the trees that it is oob on.
				observationsThatAreOob.add(i);
				ImmutableTwoValues<Double, Map<Integer, String>> predReults = predict(processedData, observationsToPredict,
						oobObservations.get(i));
				cumulativeErrors += predReults.first;
			}
		}
		this.oobErrorEstimate = cumulativeErrors / observationsThatAreOob.size();
	}

	public ImmutableTwoValues<Double, Map<Integer, String>> predict(ProcessDataForGrowing predData)
	{
		List<Integer> observationsToPredict = new ArrayList<Integer>();
		for (int i = 0; i < predData.numberObservations; i++)
		{
			observationsToPredict.add(i);
		}
		List<Integer> treesToUseForPrediction = new ArrayList<Integer>();
		for (int i = 0; i < forest.size(); i++)
		{
			treesToUseForPrediction.add(i);
		}
		return predict(predData, observationsToPredict, treesToUseForPrediction);
	}

	public ImmutableTwoValues<Double, Map<Integer, String>> predict(ProcessDataForGrowing predData, List<Integer> observationsToPredict)
	{
		List<Integer> treesToUseForPrediction = new ArrayList<Integer>();
		for (int i = 0; i < forest.size(); i++)
		{
			treesToUseForPrediction.add(i);
		}
		return predict(predData, observationsToPredict, treesToUseForPrediction);
	}

	public ImmutableTwoValues<Double, Map<Integer, String>> predict(ProcessDataForGrowing predData,
			List<Integer> observationsToPredict, List<Integer> treesToUseForPrediction)
	{

		double errorRate = 0.0;
		Map<Integer, String> observationToClassification = new HashMap<Integer, String>();

		// Set up the mapping from observation index to predictions.
		Map<Integer, List<Object>> predictions = new HashMap<Integer, List<Object>>();  // One key for each observation being predicted. The list of objects contains one entry for each tree the observation is being predicted on.
		for (int i : observationsToPredict)
		{
			predictions.put(i, new ArrayList<Object>());
		}

		// Perform the predictions for each tree.
		for (int i : treesToUseForPrediction)
		{
			Map<Integer, Object> predictedValues = forest.get(i).predict(predData, observationsToPredict);
			for (int j : observationsToPredict)
			{
				predictions.get(j).add(predictedValues.get(j));
			}
		}

		// Determine the prediction for each observation.
		if (ctrl.isClassificationUsed)
		{
			for (int i : predictions.keySet())
			{
				List<Object> predictedValues = predictions.get(i);
				Map<String, Double> predictedClasses = new HashMap<String, Double>();  // Mapping from a class to the sum of the fractions of times the observation was classed as it.
				for (Object o : predictedValues)
				{
					Map<String, Double> predValue = (Map<String, Double>) o;
					if (predictedClasses.size() == 0)
					{
						for (String j : predValue.keySet())
						{
							predictedClasses.put(j, predValue.get(j));
						}
					}
					else
					{
						for (String j : predValue.keySet())
						{
							predictedClasses.put(j, predictedClasses.get(j) + predValue.get(j));
						}
					}
				}
				// Determine the majority classification for the observation.
				String majorityClass = "";
				double largestNumberClassifications = 0.0;
				for (String key : predictedClasses.keySet())
				{
					if (predictedClasses.get(key) > largestNumberClassifications)
					{
						majorityClass = key;
						largestNumberClassifications = predictedClasses.get(key);
					}
				}
				// Record the majority classification for the observation.
				observationToClassification.put(i, majorityClass);
			}

			// Record the error rate for all observations.
			for (int i : observationToClassification.keySet())
			{
				String predictedClass = observationToClassification.get(i);
				if (predData.responseData.get(predictedClass).get(i) != 1)
				{
					// If the classification is not correct.
					errorRate += 1.0;
				}
			}
			errorRate = errorRate / observationToClassification.size();
		}
		else
		{
			// Classification is not being used.
			// PREDICT THE VALUE
			// DETERMINE ERROR (many ways to determine error (maybe abs difference is best))
			// POSSIBLY USE SQUARED DIFFERENCE
		}

		return new ImmutableTwoValues<Double, Map<Integer,String>>(errorRate, observationToClassification);
	}

//	/**
//	 * Probabilistic prediction and classification are always used here. Returns the raw probabilistic scores of the predictions. I.e. if you have 500 trees and two
//	 * classes A and B, then for each observation to be predicted it will return a hash which gives the total weight of the prediction of A
//	 * and B, the sum of the predictive weight will have to equal the number of trees.
//	 * 
//	 * @param predData
//	 * @param observationsToPredict
//	 * @param treesToUseForPrediction
//	 * @return
//	 */
//	public Map<Integer, Map<String, Double>> predictRawResults(ProcessDataForGrowing predData, List<Integer> observationsToPredict)
//	{
//		List<Integer> treesToUseForPrediction = new ArrayList<Integer>();
//		for (int i = 0; i < forest.size(); i++)
//		{
//			treesToUseForPrediction.add(i);
//		}
//		return predictRawResults(predData, observationsToPredict, treesToUseForPrediction);
//	}
//
//	public Map<Integer, Map<String, Double>> predictRawResults(ProcessDataForGrowing predData,
//			List<Integer> observationsToPredict, List<Integer> treesToUseForPrediction)
//	{
//
//		Map<Integer, Map<String, Double>> observationToClassification = new HashMap<Integer, Map<String, Double>>();
//
//		// Set up the mapping from observation index to predictions.
//		Map<Integer, List<Object>> predictions = new HashMap<Integer, List<Object>>();  // One key for each observation being predicted. The list of objects contains one entry for each tree the observation is being predicted on.
//		for (int i : observationsToPredict)
//		{
//			predictions.put(i, new ArrayList<Object>());
//		}
//
//		// Perform the predictions for each tree.
//		for (int i : treesToUseForPrediction)
//		{
//			Map<Integer, Object> predictedValues = forest.get(i).predict(predData, true, observationsToPredict);
//			for (int j : observationsToPredict)
//			{
//				predictions.get(j).add(predictedValues.get(j));
//			}
//		}
//
//		// Determine the prediction for each observation.
//		for (int i : predictions.keySet())
//		{
//			List<Object> predictedValues = predictions.get(i);
//			Map<String, Double> predictedClasses = new HashMap<String, Double>();  // Mapping from a class to the sum of the fractions of times the observation was classed as it.
//			for (Object o : predictedValues)
//			{
//				Map<String, Double> predValue = (Map<String, Double>) o;
//				if (predictedClasses.size() == 0)
//				{
//					for (String j : predValue.keySet())
//					{
//						predictedClasses.put(j, predValue.get(j));
//					}
//				}
//				else
//				{
//					for (String j : predValue.keySet())
//					{
//						predictedClasses.put(j, predictedClasses.get(j) + predValue.get(j));
//					}
//				}
//			}
//
//			// Record the majority classification for the observation.
//			observationToClassification.put(i, predictedClasses);
//		}
//
//		return observationToClassification;
//	}

}
