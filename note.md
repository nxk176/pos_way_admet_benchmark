Chạy thêm bindingdb thêm được ít thông tin
Tagonist hay antaginist = 1 or 0, ko xác định được thì =0.5

Từng protein ra được 1 loạt các smiles khác nhau, xem file smiles/master sẽ có nhiều kiể đo lường khác nhua như ic50 và kd/ki, có những cso chuyên về đo độ bám binding/activity, độ khác nhau mang ý nghĩa khác nhau, phải phân ra, chia ra thành 3 group đo khác nhau như binding/activity, với cùng protein sẽ có smiles khác nhau, ra vài giá trị khác nhua, xử lý ntn đó để những thằng smiles nào đó có tốt hơn smiles khác hay không, phải thiết kế 1 rukes để so sánh các smiles khác nhau, ví dụ có 2 con số khác nhau với cùng everrything else -> xử lý bằng average, biến thành log, trung bình? Median?

Với 1 protein này, có bao nhiêu smile, bind, stronger activity, dataset có giá trị, ví duj dùng 1 con để làm thì phải dùng 1 con xịn tương đương để xme có vấn đề gì hay không? 

1 chuối xuất phát -> 2 chuỗi: 1 chuỗi tốt hơn/1 chuỗi kém hơn, nhiều kiểu độ đo, nhiều kiểu đo lường khác nhau -> cần chuẩn hóa, xem có tương đồng với nhau hay ko, mấy cái tính chất khác thì tự add vào sau cũng đc, tự bịa ra cũng ko sao

Thường khí mà họ làm thiết kế thuốc -> đi từ cái đã biết (screening) 100k chất-> vài chất r làm thí nghiệm -> tạo thành ra 1 group lấy ra để học

Data chia ra negative/positive , ví dụ có 2 chục chất khác nhau -> 5,7 nhóm khác nhau, có tốt có xấu. positive/negative theo kiểu binding thì khó tìm đc negative, vì paper thường họ cũng ko report, hard negative cũng hiếm thấy hơn, 

Có thể lấy mô hình mô phỏng để xem negative.

Ntn là đủ dderr làm theo kiểu ranking, có cái tốt hơn, có cái xấu hơn,


Bindingdb download 1 snapshot về -9gb
Chuyển sang bài toán có engagement hay không -> bind, gì đấy, nhiều kiểu khác nhau
Khuyên nên làm sang hướng 2, vì nếu đào sâu vào hướng 1 chính là giải quyết hướng 2. Trong assay có các thông tin sâu hơn như về protein