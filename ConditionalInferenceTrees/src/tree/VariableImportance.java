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
import java.util.Set;

/**
 * @author Simon  Bull
 *
 */
public class VariableImportance
{

	/**
	 * Calculate the conditional importance of the covariables.
	 * 
	 * @param forest
	 * @param ctrl
	 * @param threshold
	 * @param numberPermutations
	 * @param calculateUsingAUC
	 * @return : A mapping from the covariables to their conditional importance.
	 */
	public Map<String, Double> conditionalVariableImportance(Forest forest, TreeGrowthControl ctrl, double threshold, int numberPermutations,
			boolean calculateUsingAUC)
	{

		if (threshold <= 0 || threshold >= 1)
		{
			System.out.format("The threshold must be between 0 and 1 exclusive. You're value is %f\n", threshold);
			System.exit(0);
		}

		ProcessDataForGrowing dataUsedForGrowing = forest.processedData;

		if (ctrl.isClassificationUsed)
		{
			if (calculateUsingAUC && dataUsedForGrowing.responseData.keySet().size() > 2)
			{
				// You can only use binary classification with AUC.
				System.out.println("AUC importance works only wih binary classification. Error rate used instead.");
				calculateUsingAUC = false;
			}
		}

		// Setup the mapping that will record the difference between the oob and permuted oob error rates
		// for each covariable for each tree.
		Map<String, double[]> permutedErrorsForTrees = new HashMap<String, double[]>();
		for (String s : dataUsedForGrowing.covariableData.keySet())
		{
			double[] tmpDoubleArray = new double[forest.forest.size() * numberPermutations];
			Arrays.fill(tmpDoubleArray, 0.0);
			permutedErrorsForTrees.put(s, tmpDoubleArray);
		}

		for (int i = 0; i < forest.forest.size(); i++)
		{
			ConditionalInferenceTree currentTree = forest.forest.get(i);

			// Determine the indices of the observations that are oob on this tree
			List<Integer> oobOnThisTree = new ArrayList<Integer>();
			for (int j = 0; j < dataUsedForGrowing.numberObservations; j++)
			{
				// For every observation, determine whether the observation is oob for this tree.
				if (forest.oobObservations.get(j).contains(i))
				{
					// If the list of trees that observation index j is oob on contains tree with index i,
					// then add j to the list of observations that are oob on this tree.
					oobOnThisTree.add(j);
				}
			}

			// Classify and determine the error rate for the observations that are oob on this tree.
			Map<Integer, Object> predictions = currentTree.predict(dataUsedForGrowing, oobOnThisTree);
			double oobErrorRate = 0.0;
			if (calculateUsingAUC)
			{
				if (ctrl.isClassificationUsed)
				{
					oobErrorRate = errorClassificationAUC(predictions, dataUsedForGrowing.responseData);
				}
				else
				{
					oobErrorRate = errorRegressionAUC(predictions, dataUsedForGrowing.responseData);
				}
			}
			else
			{
				if (ctrl.isClassificationUsed)
				{
					oobErrorRate = errorClassificationErrorRate(predictions, dataUsedForGrowing.responseData);
				}
				else
				{
					oobErrorRate = errorRegressionErrorRate(predictions, dataUsedForGrowing.responseData);
				}
			}

			// Loop through all the covariables that are split on in the tree.
			for (String s : currentTree.nonTerminalCovariablesInTree())
			{
				for (int p = 0; p < numberPermutations; p++)
				{
					// Create the list of variables to condition on.
					List<String> varsToConditionOn = createConditionalList(threshold, s, new ProcessDataForGrowing(dataUsedForGrowing));

					// Geneate a permutation of the oob observations.
					List<Integer> permutation = null;
					if (varsToConditionOn.isEmpty())
					{
						// If there are no variables to condition on, then return a permutation of the
						// indices of the oob observations.
						permutation = new ArrayList<Integer>();
						permutation.addAll(oobOnThisTree);
						Collections.shuffle(permutation);
					}
					else
					{
						// If there are variables to condition on, then determine a conditional permutation
						// of the oob observations.
						permutation = conditionalPermutation(varsToConditionOn, currentTree, dataUsedForGrowing,
								oobOnThisTree);
					}

					// Alter the oob observation covariable data records for the current covariable.
					// If original oob was [2,6,23,34] and permutation is [23,2,34,6], then
					// the value of the current covariable in observation 2 is set to the value of
					// the covariable in observation 23. Basically swap around the current covariable
					// values in the oob observations according to the permutation.
					List<Object> oldCovarValues = new ArrayList<Object>();
					oldCovarValues.addAll(dataUsedForGrowing.covariableData.get(s));
					List<Object> newCovarValues = new ArrayList<Object>();

					Object[] tmpCovarValues = new Object[oldCovarValues.size()];
					for (int j = 0; j < oobOnThisTree.size(); j++)
					{
						int oldIndex = oobOnThisTree.get(j);
						int newIndex = permutation.get(j);
						tmpCovarValues[newIndex] = oldCovarValues.get(oldIndex);
					}

					for (int j = 0; j < oldCovarValues.size(); j++)
					{
						if (tmpCovarValues[j] == null)
						{
							// If the observation is not oob.
							newCovarValues.add(oldCovarValues.get(j));
						}
						else
						{
							newCovarValues.add(tmpCovarValues[j]);
						}
					}

					// Calculate the oob error rate on the permuted observations.
					dataUsedForGrowing.covariableData.put(s, newCovarValues);  // Put the new permuted data in.
					Map<Integer, Object> permPredictions = currentTree.predict(dataUsedForGrowing, oobOnThisTree);
					double permOobErrorRate = 0.0;
					if (calculateUsingAUC)
					{
						if (ctrl.isClassificationUsed)
						{
							permOobErrorRate = errorClassificationAUC(permPredictions, dataUsedForGrowing.responseData);
						}
						else
						{
							permOobErrorRate = errorRegressionAUC(permPredictions, dataUsedForGrowing.responseData);
						}
					}
					else
					{
						if (ctrl.isClassificationUsed)
						{
							permOobErrorRate = errorClassificationErrorRate(permPredictions, dataUsedForGrowing.responseData);
						}
						else
						{
							permOobErrorRate = errorRegressionErrorRate(permPredictions, dataUsedForGrowing.responseData);
						}
					}
					permutedErrorsForTrees.get(s)[p + (i * numberPermutations)] = permOobErrorRate - oobErrorRate;
					dataUsedForGrowing.covariableData.put(s, oldCovarValues);  // Reset the data for the current covariable.
				}
			}
		}

		// Calculate the mean of the error rate differences for ecah covariable.
		Map<String, Double> returnValue = new HashMap<String, Double>();
		for (String s : permutedErrorsForTrees.keySet())
		{
			double meanValue = 0.0;
			for (double d : permutedErrorsForTrees.get(s))
			{
				meanValue += d;
			}
			meanValue /= permutedErrorsForTrees.get(s).length;
			returnValue.put(s, meanValue);
		}

		return returnValue;

	}

	/**
	 * @param varsToConditionOn
	 * @param currentTree
	 * @param dataUsedForGrowing
	 * @param observationsThatAreOob
	 * @return
	 */
	List<Integer> conditionalPermutation(List<String> varsToConditionOn, ConditionalInferenceTree currentTree,
			ProcessDataForGrowing dataUsedForGrowing, List<Integer> observationsThatAreOob)
	{

		Map<String, List<Double>> permsToPerform = new HashMap<String, List<Double>>();
		for (String s : varsToConditionOn)
		{
			if (!currentTree.nonTerminalCovariablesInTree().contains(s))
			{
				// If the covariable to condition on was not used to split a node in the current tree,
				// then continue onto the next node.
				continue;
			}

			// Determine what the split points are for the covariable.
			List<Double> cutpoints = new ArrayList<Double>();
			cutpoints.addAll(currentTree.cutpoints(s));
			if (!dataUsedForGrowing.variableTypeMapping.get(s).equals("c"))
			{
				// Don't sort the categorical cutpoints as the order is important.
				Collections.sort(cutpoints);
			}
			permsToPerform.put(s, cutpoints);
		}

		List<Integer> returnValue = new ArrayList<Integer>();
		if (permsToPerform.keySet().isEmpty())
		{
			// If there are no conditioning variables used in the tree.
			returnValue.addAll(observationsThatAreOob);
			Collections.shuffle(returnValue);
		}
		else
		{
			// If there are conditioning variables that are used to split a node in the tree.

			// Generate a list that contains one entry for every possible combination of the
			// cutpoints for the covariables that are in the tree and to be conditioned on.
			// For example, if the covariables are 'a' and 'b', and have cutpoints [1,4] and [5,7,9]
			// respectively, then there will be a list for when the values for 'a' and 'b' in the observation are:
			// 'a' <= 1 && 'b' <= 5, 'a' <= 1 && 'b' > 5 && 'b' <= 7, 'a' <= 1 && 'b' > 7 && 'b' <= 9, 'a' <= 1 && 'b' > 9
			// 'a' <= 4 && 'b' <= 5, 'a' <= 4 && 'b' > 5 && 'b' <= 7, 'a' <= 4 && 'b' > 7 && 'b' <= 9, 'a' <= 4 && 'b' > 9
			// 'a' > 4 && 'b' <= 5, 'a' > 4 && 'b' > 5 && 'b' <= 7, 'a' > 4 && 'b' > 7 && 'b' <= 9, 'a' > 4 && 'b' > 9
			// This will form a design matrix where the index of an observation is only in the list corresponding to
			// the boolean condition that is true. For example, if obs i has ('a' <= 4 && 'b' > 7 && 'b' <= 9) == true,
			// then the list for that boolean condition will contain index i.
			// The method for ordering the indices for the oobObsThatMeetCriteria lists are as follows:
			// A variable < the minimum cutpoint is given a 0, between the first and second cutpoints a 1, ..., greater than the last cutpoint
			// is given a #cutpoints + 1. E.g. 'a' < 1 is given 0, 1 <= 'a' < 4 is given 1 and a > 4 is given 3.
			// The first variable in the order is given a modifier of 1, the second variable a modifier of 2, etc.
			// Therefore, 1 <= 'a' < 4 && 7 <= 'b' < 9 would be (1 * 1) + (2 * 2) = the list with the fifth index.
			// For categorical variables, the value assigned to each category is dependent on the ordering f the categories.
			// After the categories have been recalculated, the new categories are ordered (e.g. [1,4,] 2 and 3 were mapped to 4) and
			// then if the variable has category 1 the variable is given value 0, if category 4 then 1 is given as the value.
			List<List<Integer>> oobObsThatMeetCriteria = new ArrayList<List<Integer>>();
			List<List<Integer>> oobObsThatMeetCriteriaPermuted = new ArrayList<List<Integer>>();
			int numberOfConditions = 1;
			// Determine the number of possible combinations of the cutpoints.
			for (String s : permsToPerform.keySet())
			{
				if (!dataUsedForGrowing.variableTypeMapping.get(s).equals("c"))
				{
					numberOfConditions *= permsToPerform.get(s).size() + 1; // + 1 as if you have x cutpoints, then you have x + 1 ranges that a value can fall in.
				}
			}
			for (int i = 0; i < numberOfConditions; i++)
			{
				// Add one list for each condition that could occur. Each list is used to record the observations that meet the cutpoint
				// criteria that the list represents.
				// More lists will be added later when the number of categories is known for the categorical variables.
				oobObsThatMeetCriteria.add(new ArrayList<Integer>());
				oobObsThatMeetCriteriaPermuted.add(new ArrayList<Integer>());
			}
			// Go through the oob observations, and determine which criteria each one meets.
			// Use this to determine which list to enter the observation in.
			for (int i = 0; i < observationsThatAreOob.size(); i++)
			{
				int currentObsIndex = observationsThatAreOob.get(i);
				int currentIndexOfObservation = 0;  // Used to record which criteria list the observation should be entered in.
				int numberCovarsChecked = 1;  // Used to correctly increment the currentIndexOfObservation variable.
				for (String s : permsToPerform.keySet())
				{
					// Determine which criterion the current observation satisfies.
					if (dataUsedForGrowing.variableTypeMapping.get(s).equals("c"))
					{
						// If the s covariable is a categorical covariable.
						int numberOfCategoriesForS = dataUsedForGrowing.categoricalVariableLevels.get(s);
						int numberOfCutpoints = permsToPerform.get(s).size() / numberOfCategoriesForS;
						double rowSums[] = new double[numberOfCategoriesForS];
						for (int j = 0; j < numberOfCategoriesForS; j++)
						{
							for (int k = 0; k < numberOfCutpoints; k++)
							{
								rowSums[j] += permsToPerform.get(s).get(j + (k * numberOfCategoriesForS));
							}
						}
						Map<Double, List<Integer>> duplicates = new HashMap<Double, List<Integer>>();
						for (int j = 0; j < numberOfCategoriesForS; j++)
						{
							double sumOfRows = rowSums[j];
							if (duplicates.containsKey(sumOfRows))
							{
								duplicates.get(sumOfRows).add(j + 1);
							}
							else
							{
								List<Integer> tmpList = new ArrayList<Integer>();
								tmpList.add(j + 1);
								duplicates.put(sumOfRows, tmpList);
							}
						}
						List<Integer> duplicateCategories = new ArrayList<Integer>();
						for (Double j : duplicates.keySet())
						{
							if (duplicates.get(j).size() > 1)
							{
								duplicateCategories.addAll(duplicates.get(j));
							}
						}
						Collections.sort(duplicateCategories);
						List<int[]> fuse = new ArrayList<int[]>();
						for (int j : duplicateCategories)
						{
							for (int k : duplicateCategories)
							{
								if (k < j)
								{
									continue;
								}
								else
								{
									boolean isAllTheSame = true;
									for (int l = 0; l < numberOfCutpoints; l++)
									{
										double jValue = permsToPerform.get(s).get((j - 1) + (l * numberOfCutpoints));
										double kValue = permsToPerform.get(s).get((k - 1) + (l * numberOfCutpoints));
										if (jValue != kValue)
										{
											isAllTheSame = false;
										}
									}
									if (isAllTheSame)
									{
										int tmpArray[] = new int[2];
										tmpArray[0] = j;
										tmpArray[1] = k;
										fuse.add(tmpArray);
									}
								}
							}
						}
						Map<Integer, Integer> currentCovarData = new HashMap<Integer, Integer>();
						for (int j = 0; j < dataUsedForGrowing.numberObservations; j++)
						{
							Integer observationValue = (Integer) dataUsedForGrowing.covariableData.get(s).get(j);
							currentCovarData.put(j, observationValue);
						}
						int newLevel = numberOfCategoriesForS + 1;
						List<Integer> availableFusions = new ArrayList<Integer>();
						for (int j = 0; j < fuse.size(); j++)
						{
							availableFusions.add(j);
						}
						for (int j = 0; j < numberOfCategoriesForS; j++)
						{
							if (availableFusions.isEmpty())
							{
								break;
							}
							// Get the fuse indices that have the first element of the pair they hold == j;
							Set<Integer> indicesToRecord = new HashSet<Integer>();
							for (Integer k : availableFusions)
							{
								if (fuse.get(k)[0] == j)
								{
									indicesToRecord.add(fuse.get(k)[0]);
									indicesToRecord.add(fuse.get(k)[1]);
								}
							}
							// Remove the indices found.
							availableFusions.removeAll(indicesToRecord);
							// Alter the recorded categories.
							for (int k = 0; k < dataUsedForGrowing.numberObservations; k++)
							{
								if (indicesToRecord.contains(currentCovarData.get(k)))
								{
									// If the category of the observation is to be overwritten by another.
									currentCovarData.put(k, newLevel);
								}
							}
							newLevel += 1;
						}
						Set<Integer> categoriesPresent = new HashSet<Integer>();
						for (Integer j : currentCovarData.keySet())
						{
							categoriesPresent.add(currentCovarData.get(j));
						}
						if (i == 0)
						{
							// Only add the extra lists once, not once for each observation.
							int conditionsPresent = oobObsThatMeetCriteria.size();
							for (Integer j : categoriesPresent)
							{
								// Add one list to the design matrix for each new combination that can be introduced.
								// The number of new conditions introduced is = # old conditions * number of categories.
								// If you have 3 old conditions, then a categorical variable with 2 categories brings this to 6 conditions.
								for (int k = 0; k < conditionsPresent; k++)
								{
									oobObsThatMeetCriteria.add(new ArrayList<Integer>());
									oobObsThatMeetCriteriaPermuted.add(new ArrayList<Integer>());
								}
							}
						}
						// Order the categories that are present now.
						List<Integer> cats = new ArrayList<Integer>(categoriesPresent);
						Collections.sort(cats);
						currentIndexOfObservation += numberCovarsChecked * cats.indexOf(currentCovarData.get(currentObsIndex));
					}
					else
					{
						double observationValue = (Double) dataUsedForGrowing.covariableData.get(s).get(currentObsIndex);
						double previousValue = -Double.MAX_VALUE;  // The previous value of the cutpoint. Used as the lower bound to see if an observation's value falls between two cutpoints.
						boolean cutpointFound = false;
						for (int j = 0; j < permsToPerform.get(s).size() && !cutpointFound; j++)
						{
							double currentValue = permsToPerform.get(s).get(j);  // The current (upper bound) value for the cutpoint.
							// For all the ordered recorded cutpoints.
							if (previousValue <= observationValue && currentValue > observationValue)
							{
								// If the value of the covariable for the observaton is in between the current
								// cutpoint and the last cutpoint.
								currentIndexOfObservation += numberCovarsChecked * j;
								cutpointFound = true;
							}
						}
						if (!cutpointFound)
						{
							// Then the observed value is greater than all the cutpoints.
							currentIndexOfObservation += numberCovarsChecked * permsToPerform.get(s).size();
						}
					}
					numberCovarsChecked += 1;
				}
				// Add observation i to the list representing the criterion that i fits.
				oobObsThatMeetCriteria.get(currentIndexOfObservation).add(i);
				oobObsThatMeetCriteriaPermuted.get(currentIndexOfObservation).add(i);
			}

			Integer[] tmpReturnValue = new Integer[observationsThatAreOob.size()];
			for (int i = 0; i < oobObsThatMeetCriteriaPermuted.size(); i++)
			{
				Collections.shuffle(oobObsThatMeetCriteriaPermuted.get(i));
			}
			for (int i = 0; i < oobObsThatMeetCriteriaPermuted.size(); i++)
			{
				for (int j = 0; j < oobObsThatMeetCriteriaPermuted.get(i).size(); j++)
				{
					// Move the oob observation from the original index to the new index location.
					int originalOobIndex = oobObsThatMeetCriteria.get(i).get(j);
					int originalOobValue = observationsThatAreOob.get(originalOobIndex);
					int permutedOobIndex = oobObsThatMeetCriteriaPermuted.get(i).get(j);
					tmpReturnValue[permutedOobIndex] = originalOobValue;
				}
			}
			returnValue = Arrays.asList(tmpReturnValue);
		}

		return returnValue;

	}

	/**
	 * Determines the variables that the current covariable should be conditioned on.
	 * 
	 * @param threshold
	 * @param currentCovariable
	 * @param originalData
	 * @return
	 */
	List<String> createConditionalList(double threshold, String currentCovariable, ProcessDataForGrowing originalData)
	{

		List<String> toConditionOn = new ArrayList<String>();

		// Set up the control for growing the conditioning stump.
		TreeGrowthControl conditioningCtrl = new TreeGrowthControl();
		conditioningCtrl.isStump = true;
		conditioningCtrl.isBonferronniUsed = false;
		conditioningCtrl.minCriterion = -Double.MAX_VALUE;  // Made this to be unbiased so that you can get a variable importance with any number of covariables.

		// Determine whether the new stump tree should be regression or classification.
		if (originalData.variableTypeMapping.get(currentCovariable).equals("c"))
		{
			// If the current covariable is categorical, then the stump tree will be a classification tree.
			conditioningCtrl.isClassificationUsed = true;
			originalData.responseData = new HashMap<String, List<Double>>();
			List<List<Double>> tmpResponseData = new ArrayList<List<Double>>();
			List<Integer> levelNames = new ArrayList<Integer>();
			int numberOfLevels = originalData.categoricalVariableLevels.get(currentCovariable);
			for (int i = 0; i < numberOfLevels; i++)
			{
				// One list for every level of the categorical variable.
				tmpResponseData.add(new ArrayList<Double>());
			}

			// Convert the covariable categorical information to the format for the response variable.
			for (Object o : originalData.covariableData.get(currentCovariable))
			{
				Integer currentCovarValue = (Integer) o;
				int levelArrayPosition = 0;
				if (levelNames.contains(currentCovarValue))
				{
					// If the level has already been seen before, then just get the index of the level.
					levelArrayPosition = levelNames.indexOf(currentCovarValue);
				}
				else
				{
					// If the level has never been seen before, then add the level and get the last index.
					levelNames.add(currentCovarValue);
					levelArrayPosition = levelNames.size() - 1;
				}

				// Update the temporary response data lists to reflect the current observation.
				for (int i = 0; i < numberOfLevels; i++)
				{
					if (i == levelArrayPosition)
					{
						// If the current index is that of the level of the covariable or the current observation.
						tmpResponseData.get(i).add(1.0);
					}
					else
					{
						// Not the index for the current observation's level.
						tmpResponseData.get(i).add(0.0);
					}
				}
			}

			// Update the response variable information.
			for (Integer i : levelNames)
			{
				String nameOfClass = i.toString();
				originalData.responseData.put(nameOfClass, tmpResponseData.get(levelNames.indexOf(i)));
			}
			// Update the covariable data by removing the current categorical covariable.
			originalData.covariableData.remove(currentCovariable);
		}
		else
		{
			// The current covariable is not categorical, and therefore the stump tree will be a regression tree.
			conditioningCtrl.isClassificationUsed = false;
			originalData.responseData = new HashMap<String, List<Double>>();
			originalData.responseData.put(currentCovariable, new ArrayList<Double>());
			// Make the current covariable data to be the response data.
			for (Object o : originalData.covariableData.get(currentCovariable))
			{
				originalData.responseData.get(currentCovariable).add((Double) o);
			}
			// Remove the current covariable from the list of predictor variables.
			originalData.covariableData.remove(currentCovariable);
		}
		ConditionalInferenceTree stumpTree = new ConditionalInferenceTree(originalData, conditioningCtrl);

		// Determine covariables that have a criterion value > threshold.
		for (String s : originalData.covariableData.keySet())
		{
			if (stumpTree.condInfTree.pValues.get(s) > threshold)
			{
				toConditionOn.add(s);
			}
		}
		 
		return toConditionOn;

	}

	/**
	 * @param predictions - A mapping from observation indices to their prediction.
	 * @param responseData - A mapping from the response variable to the observed values.
	 * @return
	 */
	double errorClassificationAUC(Map<Integer, Object> predictions, Map<String, List<Double>> responseData)
	{
		return 0.0;
	}

	double errorClassificationErrorRate(Map<Integer, Object> predictions, Map<String, List<Double>> responseData)
	{
		double errorRate = 0.0;  // The error rate of the predictions.

		for (int i : predictions.keySet())
		{
			// Calculate the predicted class (the majority class) for the observation.
			String maxPredClass = "";
			double maxPred = 0.0;
			Map<String, Double> classPred = ((Map<String, Double>) predictions.get(i));
			for (String s : classPred.keySet())
			{
				if (classPred.get(s) > maxPred)
				{
					maxPredClass = s;
					maxPred = classPred.get(s);
				}
			}

			// Determine the correct class for the observation.
			String realClass = "";
			for (String s : responseData.keySet())
			{
				if (responseData.get(s).get(i) == 1.0)
				{
					// This is the real class of the observation.
					realClass = s;
				}
			}

			// If the predicted and correct classes are not the same, then an error has occurred.
			if (!maxPredClass.equals(realClass))
			{
				errorRate += 1.0;
			}
		}

		// Divide the number of errors by the number of predictions.
		errorRate /= predictions.size();
		return errorRate;
	}

	double errorRegressionAUC(Map<Integer, Object> predictions, Map<String, List<Double>> responseData)
	{
		return 0.0;
	}

	double errorRegressionErrorRate(Map<Integer, Object> predictions, Map<String, List<Double>> responseData)
	{
		return 0.0;
	}

}
