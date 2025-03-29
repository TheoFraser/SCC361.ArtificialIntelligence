clear all
load dataset-letters.mat
figure(1);
images = double(dataset.images); %convert images to a double
for i = 1:12
    % This for loop iterates i from 1 to 12. In the for loop it generates a
    % random Integer from 1 to 26000 then gets a random image from dataset
    % using the random Integer. Then it resizes the image from 1x728 pixels
    % to 28x28 pixels. Then it plots the image of figure 1 with its
    % corresponding label.
    randomNumber = randi(26000);
    randomImage = images(randomNumber,:);                              
    randomLabels = dataset.labels(randomNumber);
    resized = reshape(randomImage,[28,28]);
    subplot(3,4,i); imshow(resized); title(dataset.key(randomLabels))   
end;

%This firstly sets randomPerm to randperm(26000). randperm(260000) 
%generates a 1x26000 matrix with numbers from 1 to 26000 with no 
%repeating values. Then it splits the randomPerm matrix in two and the 
%images in the dataset with the same index as randomPerm is put into
%trainingData and the labels with the same index as randomPerm is put into
%trainingLabels. The other half of randomPerm 
% is put into testData and testLabels
randomPerm = randperm(26000);
trainingData = dataset.images(randomPerm(1:13000),:);
trainingLabels = dataset.labels(randomPerm(1:13000), :);
testData = dataset.images(randomPerm(13001:26000), :);
testLabels = dataset.labels(randomPerm(13001:26000), :);

%This sets the variable trfeatures to trainingData and tefeatures to
%testData. Then it sets trlabels to a categorical array of trainingLabels.
%Lastly it sets telabels to a categorical array of testLabels
trfeatures = trainingData;
trlabels = categorical(trainingLabels);
tefeatures = testData;
telabels = categorical(testLabels);

k = 1;

[accuracyE, timerE] = euclidean(trfeatures, trlabels, tefeatures, telabels, k);
[accuracyB, timerB] = braycurtis(trfeatures, trlabels, tefeatures, telabels, k);
%[accuracyCosine, timerCosine] = cosine(trfeatures, trlabels, tefeatures, telabels, k);
[accuracyTree, timerTree] = ftree(trfeatures,trlabels,tefeatures,telabels);
[accuracyecoc, timerecoc] = fecoc(trfeatures,trlabels,tefeatures,telabels);

function [accuracy,timer] = euclidean(trfeatures, trlabels, tefeatures, telabels, k)
% Firstly it initialize an empty categorical array to store the predicted 
% labels for testing samples. Then it loops through each testing sample 
% and calculates the Euclidean distance with all training samples. After 
% that it finds the k-nearest neighbors by sorting the Euclidean distance 
% values. Then it determines the mode of the labels of the k-nearest 
% neighbors as the predicted label for the testing sample. 
% Afterwards it calculates the classification accuracy by comparing the 
% predicted labels with the true labels. Lastly it displays a confusion chart 
% to show the classification results.
    tepredict = categorical.empty(size(tefeatures,1),0);
    tic,for i = 1:size(tefeatures,1) 
        comp1 = trfeatures;
        comp2 = repmat(tefeatures(i,:), [size(trfeatures,1),1]);
        l2 = sqrt(sum((comp1-comp2).^2,2));
        [~,ind] = sort(l2);
        ind = ind(1:k);
        labs = trlabels(ind);
        tepredict(i,1) = mode(labs);
    end,timer = toc;
    figure('Name',"Euclidean distance"); 
    correct_predictions = sum(telabels==tepredict);
    accuracy = correct_predictions /size(telabels,1);
    confusionchart(telabels,tepredict);
end 

function [accuracy,timer] = cosine(trfeatures, trlabels, tefeatures, telabels, k)
% Firstly it initialize an empty categorical array to store the predicted 
% labels for testing samples. Then it loops through each testing sample 
% and calculate cosine distance with all training samples. After 
% that it finds the k-nearest neighbors by sorting the cosine distance values.
% Then it determines the mode of the labels of the k-nearest neighbors as 
% the predicted label for the testing sample. Afterwards it calculates the 
% classification accuracy by comparing the predicted labels with the true 
% labels. Lastly it displays a confusion chart to show the 
% classification results.
    tepredict = categorical.empty(size(tefeatures,1),0);
    tic,for i = 1:size(tefeatures,1)
        comp1 = trfeatures;
        comp2 = repmat(tefeatures(i,:), [size(trfeatures,1),1]);
        cos = 1 - (sum(comp1 .* comp2,2) ./ (sqrt(sum(comp1.^2,2)) .* sqrt(sum(comp2.^2,2))));
        [~,ind] = sort(cos);
        ind = ind(1:k);
        labs = trlabels(ind);
        tepredict(i,1) = mode(labs);
    end,timer = toc;
    figure('Name',"Bray Curtis"); 
    correct_predictions = sum(telabels==tepredict);
    accuracy = correct_predictions /size(telabels,1);
    confusionchart(telabels,tepredict);
end 

function [accuracy,timer] = braycurtis(trfeatures, trlabels, tefeatures, telabels, k)
% Firstly it initialize an empty categorical array to store the predicted 
% labels for testing samples. Then it loops through each testing sample 
% and calculate Bray-Curtis dissimilarity with all training samples. After 
% that it finds the k-nearest neighbors by sorting the dissimilarity values.
% Then it determines the mode of the labels of the k-nearest neighbors as 
% the predicted label for the testing sample. Afterwards it calculates the 
% classification accuracy by comparing the predicted labels with the true 
% labels. Lastly it displays a confusion chart to show the 
% classification results.
    tepredict = categorical.empty(size(tefeatures,1),0);
    tic,for i = 1:size(tefeatures,1)
        comp1 = trfeatures;
        comp2 = repmat(tefeatures(i,:), [size(trfeatures,1),1]);
        BC = 1 - 2*(sum(min(comp1,comp2),2)./(size(comp1,2) + size(comp2,2)));
        [~,ind] = sort(BC);
        ind = ind(1:k);
        labs = trlabels(ind);
        tepredict(i,1) = mode(labs);
    end,timer = toc;
    figure('Name',"Bray Curtis"); 
    correct_predictions = sum(telabels==tepredict);
    accuracy = correct_predictions /size(telabels,1);
    confusionchart(telabels,tepredict);
end 

function [accuracy,timer] = ftree(trfeatures, trlabels, tefeatures, telabels)
    % This function start a timer with tic then sets the variable knnmodel to
    % fitcknn(trfeatures, trlabels);. This trains the build in mablab
    % method with the training data and training labels. Then its sets the
    % variable predicted to predict(knnmodel, tefeatures). This build in
    % matlab function returns the predictions for what labels should
    % correspond with the training data using the knnmodel. Then it stops
    % the timer and sets the value equal to timer. Then it calculates how
    % many of the predictions are currect and then calculates the accuacy.
    % Lastly it makes a figure and creates a confusionchart with training
    % lables and the predicted labels
    tic
    knnmodel = fitctree(trfeatures,trlabels);
    predicted = predict(knnmodel, tefeatures);
    timer = toc;
    correct_predictions = sum(telabels == predicted);
    accuracy = correct_predictions /size(telabels,1);
    figure('Name',"FITCTREE"); 
    confusionchart(telabels,predicted);
end

function [accuracy,timer] = fecoc(trfeatures, trlabels, tefeatures, telabels)
    %This function start a timer with tic then sets the variable knnmodel to
    % fitcecoc(trfeatures, trlabels);. This trains the build in mablab
    % method with the training data and training labels. Then its sets the
    % variable predicted to predict(knnmodel, tefeatures). This build in
    % matlab function returns the predictions for what labels should
    % correspond with the training data using the knnmodel. Then it stops
    % the timer and sets the value equal to timer. Then it calculates how
    % many of the predictions are currect and then calculates the accuacy.
    % Lastly it makes a figure and creates a confusionchart with training
    % lables and the predicted labels
    tic
    knnmodel = fitcecoc(trfeatures,trlabels);
    predicted = predict(knnmodel, tefeatures);
    timer = toc;
    correct_predictions = sum(telabels == predicted);
    accuracy = correct_predictions /size(telabels,1);
    figure('Name',"FITCECOC"); 
    confusionchart(telabels,predicted);
end 


