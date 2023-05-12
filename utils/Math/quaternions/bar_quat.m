% Function that returns the bar product of a quaternion.
% @param[in] Q : Quaternion from which the bar product is to be calculated.
function Q_bar = bar_quat(Q)
    Q_bar = [Q(1) (-1)*Q(2:4)']';
end