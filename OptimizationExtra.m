function [cl_extra, cd_extra, cl_val, cd_val, best_Ap, best_Bp, best_Am, best_Bm] = OtimizationExtra(alpha_data, cl, thickness, camber, yc)

% 초기값 설정
ABp = [10, 10];  % Ap, Bp 초기값 설정
ABm = [10, 10];  % Am, Bm 초기값 설정
%lb = [0, 0];     % 하한값
%ub = [40, 40];   % 상한값

% 목적 함수 (최적화 대상 함수)
obj1 = @(ABp) ObjectiveFunctionPlus(alpha_data, cl, ABp, ABm, thickness, camber, yc);
obj2 = @(ABm) ObjectiveFunctionMinus(alpha_data, cl, ABp, ABm, thickness, camber, yc);

% 수동 탐색 시작(plus)
best_Ap = 1;
best_Bp = 1;
best_fval = obj1([best_Ap, best_Bp]);

for A = 1:40
    for B = 1:40
        new_fval = obj1([A, B]);
        if new_fval < best_fval
        best_fval = new_fval;
        best_Ap = A;
        best_Bp = B;
        end
    end
end
ABp_opt = [best_Ap, best_Bp];

% 수동 탐색 시작(minus)

best_Am = 1;
best_Bm = 1;
best_fval = obj2([best_Am, best_Bm]);

for A = 1:40
    for B = 1:40
        new_fval = obj2([A, B]);
        if new_fval < best_fval
        best_fval = new_fval;
        best_Am = A;
        best_Bm = B;
        end
    end
end
ABm_opt = [best_Am, best_Bm];

% 최적화된 Ap, Bp에 대한 polar 계산
[cl_extra, cd_extra, cl_val, cd_val] = ExtrapolationFn(alpha_data, cl, ABp_opt(1), ABp_opt(2), ABm_opt(1), ABm_opt(2), thickness, camber, yc, 0, 0);
end

function fval = ObjectiveFunctionPlus(alpha_data, cl, ABp, ABm, thickness, camber, yc)
    % 데이터의 마지막값들을 추출
    Aend1 = alpha_data(end);
    Aend2 = alpha_data(end-1);
    Aend3 = alpha_data(end-2);
    clend1 = cl(end);
    clend2 = cl(end-1);
    clend3 = cl(end-2);

    A = 0;

    % 기울기 변화범위 확인
    if abs(Aend1 - Aend2) < 0.2             % 기울기 변화가 0.2 미만일 때
        Aend2 = Aend1 - 0.5;
        Aend3 = Aend1 - 1;
        clend2 = cl(dsearchn(alpha_data', Aend2));
        clend3 = cl(dsearchn(alpha_data', Aend3));
    elseif abs(Aend1 - Aend2) > 1.6       % 기울기 변화가 1.6 초과할 때
        % 경고문 띄우기
    else                                    % 그 사이

    end

    % if (clend1 < clend2) && (clend2 < clend3)       %1
    %     % 정상 범위
    % elseif (clend3 < clend2) && (clend3 < clend1)   
    %     if abs(clend1 - clend3) < abs(clend1 - clend2)        %2
    %         % 정상 범위
    %     elseif abs(clend1 - clend3) > abs(clend1 - clend2)    %3
    %         if 3*abs(clend1 - clend2) < abs(clend2 - clend3)
    %             % 1,2만 계산하거나 보정치 넣기
    %             % Aend3 = 0;
    %         else
    %             % 정상 범위                                              
    %         end
    %     end
    % 
    % elseif (clend1 > clend2) && (clend2 > clend3)
    %     if (abs(clend3 - clend2))/2 < abs(clend2 - clend1)
    %         %그냥
    %     else
    %         %직선 경고문 띄우기
    %     end
    % end

    % alpha값을 이용하여 평판유동의 cl값 계산
    [~,~,cl_val1,~] = ExtrapolationFn(alpha_data, cl, ABp(1), ABp(2), ABm(1), ABm(2), thickness, camber, yc, Aend1, A);
    [~,~,cl_val2,~] = ExtrapolationFn(alpha_data, cl, ABp(1), ABp(2), ABm(1), ABm(2), thickness, camber, yc, Aend2, A);
    [~,~,cl_val3,~] = ExtrapolationFn(alpha_data, cl, ABp(1), ABp(2), ABm(1), ABm(2), thickness, camber, yc, Aend3, A);

    
    % 최적화 대상은 Cl의 오류 최소화
    clerror1 = (abs(clend1-cl_val1))^2;
    clerror2 = (abs(clend2-cl_val2))^2;
    clerror3 = (abs(clend3-cl_val3))^2;

    % 목적 함수는 Cl값의 합 차이로 정의
    fval = (clerror1+clerror2+clerror3)/3;
end 

function fval = ObjectiveFunctionMinus(alpha_data, cl, ABp, ABm, thickness, camber, yc)
    % 데이터의 처음값들을 추출
    A1 = alpha_data(1);
    A2 = alpha_data(2);
    A3 = alpha_data(3);
    cl1 = cl(1);
    cl2 = cl(2);
    cl3 = cl(3);

    Aend = 0;

    % 기울기 변화범위 확인
    if abs(A1 - A2) < 0.2             % 기울기 변화가 0.2 미만일 때
        A2 = A1 + 0.5;
        A3 = A1 + 1;
        cl2 = cl(dsearchn(alpha_data', A2));
        cl3 = cl(dsearchn(alpha_data', A3));
    elseif abs(A1 - A2) > 1.6       % 기울기 변화가 1.6 초과할 때
        % 경고문 띄우기
    else                                    % 그 사이

    end

    if (cl1 > cl2) && (cl2 > cl3)       %1
        % 정상 범위
    elseif (cl3 > cl2) && (cl3 > cl1)   
        if abs(cl1 - cl3) < abs(cl1 - cl2)        %2
            % 정상 범위
        elseif abs(cl1 - cl3) > abs(cl1 - cl2)    %3
            if 3*abs(cl1 - cl2) < abs(cl2 - cl3)
                % 1,2만 계산하거나 보정치 넣기
                A3 = 0;
            else
                % 정상 범위                                              
            end
        end

    elseif (cl1 < cl2) && (cl2 < cl3)
        if (abs(cl3 - cl2))/2 < abs(cl2 - cl1)
            %그냥
        else
            %직선 경고문 띄우기
        end
    end

    % alpha값을 이용하여 평판유동의 cl값 계산
    [~,~,cl_val1,~] = ExtrapolationFn(alpha_data, cl, ABp(1), ABp(2), ABm(1), ABm(2), thickness, camber, yc, Aend, A1);
    [~,~,cl_val2,~] = ExtrapolationFn(alpha_data, cl, ABp(1), ABp(2), ABm(1), ABm(2), thickness, camber, yc, Aend, A2);
    [~,~,cl_val3,~] = ExtrapolationFn(alpha_data, cl, ABp(1), ABp(2), ABm(1), ABm(2), thickness, camber, yc, Aend, A3);

    % 최적화 대상은 Cl의 오류 최소화
    clerror1 = abs(cl1-cl_val1);
    clerror2 = abs(cl2-cl_val2);
    clerror3 = abs(cl3-cl_val3);

    % 목적 함수는 Cl값의 합 차이로 정의
    fval = clerror1^2+clerror2^2+clerror3^2;
end
