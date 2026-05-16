[00:00] nó đã tổ hợp data thì từ đội thì ngay từ đầu bọ em đã thu thập từ Popcam từ
[00:07] CBL à trong tuần này thì bọn em đã làm được thêm được là Papirus với lại phải tìm được
[00:13] là có Bing DB ấy thì tổng cộng là đã tìm được khoảng 300.000 molecule mà bọn em đang tổng
[00:21] hợp lại để xem là các điều kiện nó có giống nhau không và trong thời gian đấy thì trong
[00:26] tuần vừa rồi thì em cũng đã đọc về hai bài thứ nhất là bài là hai bài về benchmark
[00:31] thì bài tượng đầu tiên nó tên là Mall L bench à được accept ở IC 2026 bài thứ hai
[00:37] cũng là một bài IC 2026 là bài MC IQ
[00:43] thì phần vừa rồi thì chủ yếu em cũng đọc về cái khoản Rattle của các bài ấ nên thế
[00:48] là có thể là có một số sai sót thì có gì em xin phép ạ
[00:58] thì đầu tiên là về bài Morum Bench thì nó là một bài một bài benchmark để đánh giá những
[01:05] các cái mô hình về thứ nhất là structure recognition thứ hai là về khoản editing và thứ ba là
[01:11] về generation thì đây là một cái pipeline mà những cái nhà hóa học mà làm đầu tiên là họ
[01:20] sẽ analyze được cái structure của cái molecule gì về scaffold về các nhóm trúc và về hình thái 3D
[01:28] và sau đó họ sẽ reason về những cái modification để có thể đạt được cái target và cuối cùng
[01:34] là sẽ apply những cái change để có thể có được cái molecule mới thì cái khi mà chúng ta
[01:42] áp dụng LM vào thì nó có ba cái vấn đề thứ nhất là black box như thế là à
[01:48] chưa chắc là chúng ta đã là đã reason đó có thể học thuộc và thứ hai là những cái
[01:53] description nó có thể chưa được rõ rõ ràng như thế là có thể theo tác giảose thì nó có
[01:59] thể tạo ra một cái vấn đề là many to one maping problem cho những cái molecule và thứ ba
[02:05] là cũng chưa có nhiều benchmark mà hiệu quả thì monum bench thì nó sẽ address được cái g này
[02:13] ạ thì trong cái benchmark này người ta sẽ có ba t khác nhau đầu tiên thì nó sẽ là
[02:21] recognition thì chúng ta có một molecule M và có một cái TQ và chúng ta sẽ à tạo ra
[02:28] được ra được một cái câu trả lời từ cái quy đó à thứ hai thì nó sẽ là molecule
[02:36] editing thì có một cái molecule thêm edit vào thì chúng ta sẽ có một cái molecule đây là molecule
[02:43] editing bình thường và thứ ba là generation thì có một cái description thì chúng ta sẽ tạo ra được
[02:48] một cái machine molecule
[02:56] thì về recognition thì chủ yếu tác giả nói về topology connectivity các nhóm chức các cái substructure và về
[03:04] hình thái 3D nó có các cái có các cái properties như sau
[03:14] còn đây là những yêu cầu liên quan đến molecule editing và molecule generation nhưng mà bình thường khi mà
[03:19] đánh giá thì tác giả lại dùng Ty Moto Similarity
[03:25] thì dataset của nó sẽ bao gồm đầu tiên là sẽ có một cái ảnh của reference molecule tiếp theo
[03:33] là sẽ có answer và sẽ có instruction ở trong máy
[03:40] thì đây là cấu trúc của dataset mà mà tác giả đã đưa trên hacking face sẽ có đầu tiên
[03:47] là sẽ có original smile sẽ có edit smile có insertion có có canary và có ảnh của original và
[03:56] edit image
[04:04] thì để đánh giá thì nó thì tác giả à propose một cái pipeline để unố những cái molecule đó
[04:13] thì đầu tiên là họ sẽ đưa ra các cái editing instruction và structural description và người ta sẽ đảm
[04:20] bảo rằng là chỉ có một mod chỉ có một modified structure thôi và thứ hai sẽ có một người
[04:27] để kiểm annoter thứ hai sẽ đánh giá lại cái anod thứ nhất và họ sẽ lặp lại cho đến
[04:36] khi nào mà hai người đều đồng ý và step 3 thì sẽ có hai người valid và sẽ check
[04:41] lại là họ sẽ không có cái ground molecule nhưng mà họ sẽ cố gắng để construct lại được cái
[04:48] structure bằng cái software ạ thì nếu mà fail thì họ sẽ hoặc là nhìn lại hoặc là bỏ cái
[04:55] hoặc là bỏ cái molecule đó đi ạ đây là minh họa của Pippine mà tác giả đã propose ở
[05:02] trong bài báo gồm ba bước ạ một là editing hai là cái review và ba là valid lại bằng
[05:09] hai người khác ạ thì ở đây tác giả bài bài benchmark thấy tác giả đã đánh giá trên bốn
[05:17] loại LMS khác nhau đầu tiên là nonreasoning tiếp theo là reasoning models và có những model vision language và
[05:25] cuối cùng là có các cái chemistry specific model thì họ đánh giá trên ba tas thì task thì họ
[05:33] sẽ dùng zero shorting và với task đầu tiên là recognition thì họ sẽ sử dụng pass one và for
[05:42] và cho cái tas là editing và generation thì họ sẽ dùng tanimoto similarity ạ ở trong ở trong repal
[05:50] thì cũng có một đoạn mà họ đã nói là tại sao à Car đã hỏi là tại sao các
[05:54] ông lại dùng smile mà không dùng selfies thì tác giả trong đây thì ở trong trong này tác giả
[06:01] đã propose và thử nghiệm trên GPT O3 nhưng mà kết quả cũng sẽ chưa cũng chưa được tốt như
[06:08] thế là tác giả cũng reasoning bảo là đây là một lỗi có thể là do LMC rõ được cái
[06:12] cấu trúc mới của selfies thì các cái key finding thì GP5 là perform tốt nhất trong cả ba t
[06:23] là về đầu tiên là cả recognition thứ hai là cả editing và generation thì đều form tốt nhất ạ
[06:34] còn đây là bảng kết quả của toàn bộ bài benchmark thì đây là đây là cái benchmark của cái
[06:44] structure recognition task được đánh giá bởi F được đánh giá bằng fast one hơi hơi hơi dài có tất
[06:52] cả các mô hình ở đây GPT4O GVT 4.5 đây này là của structure recognition có ba bản 1 2
[07:06] và trong tất cả các bài thì GP5 là được có phong tốt nhất tiếp theo thì ở bảng tiếp
[07:12] theo thì chúng ta có graph của à có table của molecule editing và generation thì ở đây trong hầu
[07:20] hết tất cả các thì5 vẫn là phong tốt nhất và sau đó là gamin 2.5 nó bài thứ nhất
[07:34] em cũng có tổng hợp một vài cái chủ yếu là em đọc về quy trình rebattle của tác giả
[07:39] và các reviewer như thế là có thể là em sẽ làm hơi hơi nhanh thầy có câu hỏi gì
[07:45] không ạ à không mấy cái magic này là cái gì đấy chưa kịp nhìn đây ạ nó sẽ có
[07:54] tas nó sẽ đây có có các t này ra thì đây là nó sẽ chia làm ba task chính
[08:01] và đây sẽ là trong này sẽ có cái subtask thì đây là task thứ nhất đây là bảng của
[08:06] task thứ nhất nó rất là dài là task recognition có các sát tác là one hop two hop thonary
[08:14] ring bond connection và hết và đây là các model khác nhau đánh giá trên cái model đó
[08:24] đ Vâng ạ đây nó giải tận ba bảng
[08:36] còn đây là hai bảng tiếp theo thì nó tính bằng Tani Motor Similarity cho hai T editing và Generation
[08:57] thế cái dataset này nó có khác gì cái dat set mà bọ em đang định làm không
[09:06] em thấy nó em cũng chưa biết được là nó có cái gì ờ thế bây giờ đây là xem
[09:16] bọn em định làm thì nó có gì hơn cái này thì bọn em đang tính là thêm những cái
[09:21] properties của molecule thì những cái tính chất lý hóa đã có sẵn rồi này những cái đấy thì cũng
[09:28] có khá nhiều rồi nhưng vậy em định thêm các cái à chỉ số nữ của weblab nữa như hôm
[09:37] nọ thầy bảo vậy ạ nghĩa là ở đây nó có instruction và target smile đúng không vâng nó có
[09:48] cả ảnh thì bọn em bọn em không có ảnh đúng không bọn em cũng là ba cái cột đầu
[09:57] ảnh ảnh thì từ SM sinh ra được rảnh à vâng nếu mà nhưng mà bọn em có ba cái
[10:02] cột đầu thế thì bọn em có cũng có ba cột đầu giống nó đúng không thế bọn em hơn
[10:07] chúng nó là cái instruction kia bọn em có thêm các cái tính chất hóa lý à mấy cái biding
[10:13] à
[10:16] mấy cái bing anh nó sẽ vâng ạ nó nó sẽ có thêm các cột em mong là nó sẽ
[10:21] thêm các cột à có thêm cột đấy ạ ngoài cái cột instruction kia thì còn có cả thêm cột
[10:28] về kiểu loại binding các thứ à hả ok th ra này dataset này thì chủ yếu nó chỉ liên
[10:38] quan đến việc sửa đổi công thức thôi đúng không vâng còn còn cái cái cột cột thứ tư này
[10:44] ý nghĩa là gì cột thứ tư này em nó em em cũng chưa tìm hiểu nó là gì thầy
[10:49] ạ nó có một kiểu một cái mã ID của nó chắc là nó sẽ mapping đến một vài cái
[10:53] data khác và em sẽ nghiên cứu thêm paper là nó không chỉ cả cái này là cái này là
[10:59] chỉ instruction là chỉ để sửa cái công thức từ một chất sang chất khác
[11:10] nhưng mà cái chất edited smile này là người ta lấy từ đâu vậy người ta lấy từ phòng thử
[11:16] nghiệm hóa hay là không hình như là nó sinh ra xem có người review xem có đúng hay không
[11:22] à vâng một chất là có bốn người review ấy ạ thế thì cũng là valuable rồi em đọc là
[11:30] lúc đầu người ta trong rep người ta hỏi là tại sao các ông làm ít chất thế thì họ
[11:36] lấy số liệu người ta làm thực tế ra chúng mình gọi là lấy những cái người ta tính được
[11:42] là làm bao nhiêu giờ và làm bao nhiêu chất như này như kia và người ta so với một
[11:47] thằng một thằng mà thằng GPT dùng để benchmark mà là có 198 thì còn ông này thì có 200
[11:55] chất lá gốc cộng với 200 cái nữa là khoảng tầm 400 cái em em quay lại em có cái
[12:03] slide nào có cái thống kê của cái datet này không bao nhiêu một chất có thể đó đây đây
[12:11] nó là các instant kiểu nó là các một cái
[12:17] em cũng chưa hình dung được rõ ràng thầy em muốn hỏi thầy nó có kiểu em cũng chưa hiểu
[12:22] là nó có cái kiểu các nguyên tử hay là cái chất for gì đấy mà thế có download được
[12:28] cái này về không có thầy ạ nó có trên hing face download về thì nó có khoảng bao nhiêu
[12:34] dòng Đây em xem nó đây trong edit nó sẽ có 400 row trong test thì nó chỉ có 200
[12:42] row à ừ thế thì ít thật
[12:52] recognition thì người ta augment lên 33,8.000 nhưng mà bởi vì nó có rất nhiều cái subas nó bé trong
[12:57] đấy tức là cái gốc cái mà có người review thật nó chỉ cóả nó chỉ có khoảng 400 chất
[13:05] đúng không 400 cấu trúc đúng không vâng ừ thế thì
[13:14] sử đâu ờ I 2026
[13:23] thuê mấy ông giáo tự hóa có clear
[13:29] năm nay năm nay thì em có tìm về cái tab liên quan đến molecule thì có bốn bài đấy
[13:36] thì ở đây cùng bày hai bài
[13:42] 400 chất thì hết tầm bao tiền hả thầy sa biết được thì
[13:55] cái này tắt 1 phút một chất chứ nhỉ em em thấy cái người ta thuê GPU để nhiều tiền
[14:03] rồi tạo data thì cần gì
[14:13] chắc là test mấy cái ô một tại vì thực ra làm cái này thì chả ai thao tác trực
[14:20] tiếp đến cái chỗ smile cả nó vẽ ra xong rồi khi mà có cái instruction thay gốc tì bằng
[14:27] gốc tì nó vẽ luôn cho cái đấy xong nó sẽ chuyển lại smart bằng phần mềm
[14:36] đây rồi thấy em có bài phát không có một bài nữa thì em cũng thấy nó ít mà họ
[14:42] defend được cái kết quả của họ như thế là em cũng thấy khá là hay thì nó dữ liệu
[14:47] thật thì đi cái đâu đấy
[15:05] ở đây có một vài thông tin về tại sao mà structure recognition thì họ chia ra làm hai tas
[15:12] thầy cũng muốn đọc thêm em có một đây họ kia làm hai tas thì đây là thứ nhất là
[15:20] recognition và thứ hai là localization thì họ lý giải tại sao là localization nó cái performance của localization nó
[15:28] lại giảm đi thì do dùng by encoding thì ví dụ những cái mà index nó gần nhau thì họ
[15:34] lại merge vào một cái token như thế là ví dụ như là item wise tokenizer thì sẽ có performance
[15:40] tốt hơn thì ở đây cũng 100 ờ nói chung là cái bài này là nó làm thêm 400 chất
[15:49] kia mới xong rồi chủ yếu data là nó lấy ở các cái bộ đã có rồi đúng không xong
[15:55] rồi xong rồi nó chạy mấy cái mấy cái downstream task như kiểu recognization rồi thì editing xong rồi generation
[16:06] gì đấy xong rồi họ thử với cảm với cái lm nữa đúng không vâng chủ yếu là họ benchmarkên
[16:11] rất trên rất nhiều con và xong rồi có nghĩa là một có một đống các cái bảng kết quả
[16:16] chứ gì vâng thế bài này điểm cao không bài này bài này em nhớ bài này là borderl borderl
[16:24] res ba ông review cho ba con bốn nhưng mà họ defend được lên hết sáu mà vẫn except ok
[16:36] rồi bài không có gì chỉ có tốn tiền tốn thời gian một chút tốn tiền thú
[16:46] em thấy nhiều kết quả quá thì nó gọi ra xong rồi nó tự phân tích kết quả thì nó
[16:53] có phải đập trình quá nhiều đâu
[16:58] và test chậm nhanh dữ liệu chỉ có 400 chất hình như nó test cả về những cái dữ liệu
[17:07] cũ nữa rồi thầy dữ liệu nó lấy từ các cái bộ khác à bộ khác thì sao nó có
[17:14] nhãn mà nó so sánh được ừ đấy Tùng
[17:21] dạ thầy đang hỏi là nếu mà bây giờ data của họ là có 400 chất mới đúng không khoảng
[17:27] tầm đấy thế còn dữ liệu còn lại là họ lấy từ cái bộ khác đúng không dữ liệu nào
[17:33] ạ tức là họ có test với dữ liệu ở ngoài 400 chất này không theo mình hiểu là không
[17:42] em chưa ấy ô thế thì cái bộ dữ liệu này chỉ có 400 chất thôi hả Tùng sau nãy
[17:47] nhìn thấy có mấy mấy triệu gì cơ mà xanh lắm nó thu mấy triệu đấy nó thu lại thì
[17:53] nó sẽ chọn ra được 400 chất thôi nó ở đoạn đầu anh ạ nó ở đoạn đầu
[18:02] đấy 200 đây này thấy 178 molecule đấy là nó lấy từ cái unik ra cái này em à ok
[18:11] 200 ok
[18:16] ừ thế ừ dữ liệu có mấy trăm chất này mà cũng public được họ họ đ chủ yếu là
[18:24] em thấy là họ defend rất đẳng cấp anh ạ chủ yếu thấy rebattle họ defend đẳng cấp thật kiểu
[18:29] bảo họ cũng so sánh với lại những cái thằng bas mà top bây giờ thì như vừa nãy họ
[18:35] họ họ chứng minh là có một thằng có một ông đây em ừ thế nhưng mà cái các cái
[18:41] những cái data set hiện tại thì nó có tốt hơn cái của của của bọn này không không ạ
[18:50] đây tốt nhất đây em nhớ nhá nó là một cái nó tên Vâng ạ cái mà đây nó tên
[18:54] là GPQA Benchmark thì nó chỉ có 198 thô để cho để dùng cho thằng JPT ok ok của cô
[19:02] sẽ 200 ông này nhiều tiền hơn hết một bài
[19:21] tiếp theo thì em sẽ present bàn molecular IQ ờ cho mình hỏi cho mình hỏi lại cái bài trước
[19:27] mỗi không không để ý kỹ lắm tức là cái prom thì họ sinh ra như nào họ từ chất
[19:31] họ sinh ra prom hay là họ từ chất cộng prom họ sinh ra chất mới ừ trong cái này
[19:38] thì em thấy em cũng chưa rõ ràng như em em đang em chỉ nghĩ là các đây họ sẽ
[19:44] là làm manual hết ờ nhưng cái instruction ấy họ sinh ra instruction như thế nào cái đấy em chưa
[19:53] nghiên cứu đến thầy ạ tại vì nếu mà nó dễ quá chỉ có thay gốc này bằng gốc kia
[19:59] thì tức là nhiều khi này nó chỉ bias vào một cái nhóm tas nào đó khá
[20:11] đây đây thấy ạ nó có một cái đây ví dụ đây ạ đây là một instruction ừ
[20:23] nhưng mà ví dụ phần đấy là thay này remove này hoặc là chặt ti này có một số cái
[20:29] thì nó ấy vâng đấy tức là làm thế nào ra được các kiểu instruction này
[20:39] tôi được biết thế à vâng
[20:57] Thì bài tiếp theo em sẽ present là đó là bài molecular IQ thì chủ yếu nó sẽ benchmark cái
[21:04] reasoning capabilities của các con với nhau
[21:12] thì có những cũng đã có những cái benchmark khác nhưng mà molecular IQ thì sẽ benchmark những cái task
[21:20] mà cái correctness list của nó có thể verify lại bằng những cái tool như là a à như như
[21:25] là adicit họ sẽ đưa ra một cái phương pháp đánh giá để có để nó sẽ gọi là chính
[21:32] xác hơn hoặc là so với những cái con khác ví dụ những cái mô hình khác nó sẽ dùng
[21:39] ví dụ có có các benchmark sẽ dùng các cái con mô hình khác để nó đánh giá lại nhưng
[21:43] mà con này thì nó sẽ thiết kế ra một cái benchmark làm sao phải kiểu nó ừ để kiểu
[21:51] mình có thể đánh giá lại kiểu định lượng ấy adic kit thì motivation của bài toán thì đầu tiên
[21:58] họ sẽ họ nói ra là ví dụ những con thì có thể học thuộc thay vì reasoning và sẽ
[22:06] không thể phân biệt được những cái pattern từ cái graph và Cái cái answerol thì có thể bị contaminate
[22:15] data
[22:18] hoặc là những cái benchmark khác thì có thể ẫn đơn nó vẫn đơn giản quá hoặc là vẫn chưa
[22:24] thể phân biệt được cái structure từ những cái label mà đã được nhớ sẵn hoặc là có thể những
[22:33] cái verifier có những cái bias khi mà đánh giá
[22:40] về khoản thứ ba về motivation thì có những cái benchmark khác thì có những cái benchmark là đã bị
[22:46] saturate tức là bị bão hòa rồi hoặc là có những cái benchmark bị contaminate
[22:57] thì về molecular IQ thì họ sẽ đánh giá trên ba tas chính thì đầu tiên là feature counting ví
[23:03] dụ như là có bao nhiêu x ở trong molecul m này có bao nhiêu bao nhiêu rings có bao
[23:07] nhiêu substructure thứ hai họ sẽ index ra những cái cái attribute thì ụ ở đây là có một câu
[23:16] hỏi at which item index is expresent đấy và sau đó sẽ có Q&A nà đó và thứ ba thì
[23:26] nó sẽ là generation thì là generate một cái smile thỏa mãn một cái constraint đã cho trước
[23:52] thì ở đây họ sẽ đánh giá ở trên có ba có ba accis để đánh giá thì đầu tiên
[23:56] sẽ đánh giá bằng smile representation thứ hai họ sẽ đánh giá trên cái complexity và thứ ba họ sẽ
[24:02] đánh giá trên một cái gọi là multitask load với n1235 thì chủ yếu hôm nay em sẽ present về
[24:09] khoản một trước khoản hai và khoản ba thì sẽ có hai bảng khác nhưng em cũng chưa hiểu về
[24:14] cái data lắm nên thế là em sẽ về cái trước thì về cái dataset construction thì đầu tiên họ
[24:22] sẽ lấy từ upcam và họ sẽ đảm bảo làm sao cho nó đều single fragment và đều có cbon
[24:29] ở trong đó thì sẽ có khoảng vài triệu candidate ờ sau đó họ sẽ đến bước filtering giữ lại
[24:36] khoảng từ 5 đến 50 heavy items và giữ cho cái length của cái chỗ smile nhỏ hơn 100 và
[24:43] làm sao cho cái kid mà nó pass được cái smiles đó tiếp theo là đến bước decontamination thì họ
[24:50] sẽ xóa họ sẽ bỏ những cái molecule mà có ở trong trong những cái evaluation set như ở bên
[24:58] trên tiếp theo họ sẽ cluster bằng cái similarity của cái structure của nó và họ sẽ chia thành các
[25:06] cái tập trang tập test tập trending pool là có 1,3 triệu molecules sẽ chia làm hai khoảng test là
[25:13] easy test và hard test sẽ mỗi cái có 1 triệu và cuối cùng họ sẽ sampling nó ra họ
[25:22] sẽ sampling 849 molecule ra và sẽ tính BS này là cái chỉ số về complexity và tiếp theo là
[25:31] họ sẽ chia thành các bin và để chia thành các cái multitask load thì tổng cộng là họ đang
[25:37] propose khoảng 51 question
[25:47] thì họ sẽ họ sẽ sử dụng các cái feature ở đây có tổng cộng có 30 có 30 feature
[25:52] trên thành sáu nhóm
[25:58] À về thứ nhất là về nhóm chức thứ hai về composition về về topology thì có cam type và
[26:03] graph type về chemical perception và cuối cùng là synthesis và fragmentation
[26:11] thì cách mà tính điểm của họ thì đầu tiên họ sẽ cho input model thì sẽ có một cái
[26:17] cũng như các model khác là họ sẽ dùng một một smile scope cộng với một task và họ sẽ
[26:23] assign smile tiếp theo là họ sẽ đánh giá bằng họ sẽ gọi các mô hình để vào để dùng
[26:33] API sẽ preprocess và extract ra được những cái đặc điểm của nó và họ sẽ tạo ra cái đáp
[26:39] án thì ở trong Reportal thì họ sẽ có các cái roll out tức là kiểu các cái option khác
[26:45] nhau thì ở đây tác giả lúc đầu trong bài trong bài báo thì propose ba roll out thế mà
[26:52] khi mà RER hỏi thì người ta lại hỏi là các ông bây giờ các ông tăng số run out
[26:58] lên thì có được không thì người ta chứng minh là tăng số run out lên thì ranking cũng đã
[27:02] không thay đổi nó vẫn reasoning vẫn làm được tốt mô hình vẫn làm được tốt nên thấy là không
[27:08] có vấn đề gì xảy ra ạ và cuối cùng là họ và tiếp theo thì họ sẽ extract ra
[27:12] được những cái thông tin à thông tin về hai à thông tin về kiểu xếp tầng thì họ sẽ
[27:19] xếp và sắp xếp các cái thông tin về molecule lại với nhau và cuối cùng là họ sẽ sử
[27:24] dụng kit để làm được cái bước scoring và và verify được đáp án thì với ở trong proposal ban
[27:36] đầu thì ví dụ trên 2/3 roll out mà ok thì họ sẽ cho đó là một success còn lại
[27:41] không thì chủ yếu trong bài này thì họ evaluate và họ đưa ra cái kết quả dựa trên accuracy
[27:50] và họ đánh giá trên khá nhiều model thì họ đánh giá trên hai trên hai nhóm chính đầu tiên
[27:55] là chemistry LM và Generalist Lom thì hầu hết các chemistry LM thì đều perform khá là tệ còn à
[28:02] các mô hình mà nhiều tham số thì tôi phong tốt hơn
[28:13] nên về độ chính xác thì à họ chia accuracy dựa trên ba substance bé là cting generation đó là
[28:21] hết của bài này rồi
[28:30] ừ ừ ok mình chưa có câu hỏi gì chưa nghĩ ra được bài này cái dataset của nó họ
[28:39] up họ đang trên git mà họ dùng chpkl em cũng chưa tích được để xem là dataset nó trông
[28:44] như nào và em sẽ xử lý thêm
[29:06] ok còn gì nữa không về về surve thì em cũng xin hết ạ
[29:29] mình phải bit review một số paper về dataset và B chưa có thời gian đọc có danh sách các
[29:37] cái title object thì có gì t mình gửi thôi nhưng mà đừng có gửi cho ai nhá xong xem
[29:42] qua xem có thấy cái idea nó hay không thì nốt lại cho mình để mình biết cái bài đấy
[29:52] Ờ cái table này là các cái table của new hôm nay
[30:02] gửi xong zalo
[30:13] đây cũng là hai cái trang đầu tiên này có nhiều phết nên cũng không chưa có thời gian print
[30:19] mấy cái trang sau mà
[30:42] đế cái trang thứ ba này nó ít về chất rồi
[30:54] đây là người ta đang ở đâu ạ thầy đây là các bài New reviewer mà của 2025 ạ năm
[31:05] nay đang chuẩn bị review thầy đang thầy đang là review của
[31:15] em chơi public lên trên overview bọn em không không gửi đi đâu nha đọc cái đọc đọc title và
[31:22] up thôi để xem vâng vâng ok có cái gì hay không mà bài nào quan tâm thì để mình
[31:28] biết cái bài đấy tên là bit thôi chứ có được gán hay không thì biết bit là cái gì
[31:34] hả thầy bit là đăng ký để review vào đấy
[31:43] còn họ có chọn hay không thì chưa biết nh quy quên không biết thế là họ gán cho những
[31:49] cái bài lâu lâu
[31:56] thế bây giờ ví dụ em xem ở trên này là sẽ có đúng không nhỉ ờ sang file PDF
[32:01] đấy thì chắc thì em chỉ xem được title cả thôi em có có nhấn vào đường link thì mình
[32:06] vào xem được thì cũng chỉ có title và object thôi chứ cũng chả có gì khác em chắc là
[32:09] người ngoài ấn vào link cũng không xem được đâu dạ ờ mà có xem thì chị được cái thông
[32:13] tin này khi mà sau đó thành reviewer rồi thì mới có thể down lại thì đây chỉ biết được
[32:23] title với object thôi đây trang thứ ba đây nhưng mà trang thứ ba thì nó ít liên quan chắc
[32:28] ba chà vâng
[32:41] nhưng mà ngoài idea làm data set ra để có cái gì muốn làm không
[32:53] vẫn em thì cũng ngoài ngoài ngoài ra
[32:58] nếu mà nếu mà ngoài đ thì như thầy cô với lại long có trao đổi thì có một hướng
[33:06] làm nữa ạ ừ thòng thuốc này có một số bài toán khá là gọi là kinh điển tức là
[33:19] ừ thế nào nhỉ
[33:24] thứ nhất là zero short hoặc là short learning ờ ví dụ từ Bing DB hoặc là CBL chẳng hạn
[33:33] thì có rất nhiều có nhiều target thì họ có thể có thí nghiệm với các chất khác nhau nhưng
[33:39] mà thí nghiệm thì nó lại có nhiều kiểu thí nghiệm khác nhau kiểu đo cái bing có kiểu đo
[33:45] cái nồng độ ức chế em sẽ thấy là nhiều cái unit khác nhau chỗ IC50 chẳng hạn hoặc là
[33:52] KD KI các thứ thì là có nhiều target và có nhiều chất có kết quả thí nghiệm nhưng mà
[34:00] khi mà thiết kế thuốc ấy thì thường là cái target đấy là mới các cái đấy thì thường là
[34:05] chưa có thí nghiệm nào cả chưa có cái thông tin cái chất nào ờ nào cả cho nên liệu
[34:13] mình có thể làm cái model để predict một chất có bài vào cái target đấy không dựa trên những
[34:19] cái chất khác có thể là những chất mà gần giống trên các chất và các cái target khác đã
[34:27] có data thế đấy là bài toán nói chung lúc nào cũng cần
[34:36] thì target thì nó lại có nhiều kiểu ví dụ như nó chỉ là tên của protein thôi nhưng từ
[34:41] tên của protein thì truy cập cái INDPR thì sẽ viết cái chỗ amino axit hoặc là truy cập Pam
[34:50] thì có thể nó sẽ có cái 3D cái Xray nữa không có 3D để dùng F thì cũng có
[34:56] thể predic được mà nó cũng có thể không chính xác thì đấy là từ chuỗi amino axit thì có
[35:02] thể ra được hình dạng 3D nhưng mà hình dạng 3D đấy nó lại có thể có nhiều cái bing
[35:07] pocket khác nhau những cái chỗ mà chất nó có thể bám vào cái protein đấy thế còn có những
[35:13] cái hốc chẳng hạn thì nó là cái bing chủ yếu nhưng mà có thể là những cái chất nó
[35:18] lại bay vào chỗ khác và khi nó bay ở chỗ khác ấ nó thay vì nó bật hoặc là
[35:23] tắt cái chức năng của protein thì nó lại có thể điều chỉnh ở cái con số giữa 0 và
[35:29] 1 dụ nó làm yếu đi chứ không phải là nó tắt hẳn thì nhiều khi người ta cũng muốn
[35:35] thiết kế thuốc theo kiểu đấy thay vì bật tet tắt người ta có thể điều chỉnh được thì những
[35:41] cái binding sai mà nó không phải cái chính kia xác định được nó hay không cũng là kiểu
[35:52] à cả
[35:55] buy chẳng hạn thì nhiều khi paper thì họ chỉ report những cái chất có bing chứ những cái chất
[36:00] mà không bao giờ bing thì không có viing thì họ không report thì nó không public được như là
[36:06] thường nếu mà lấy data trên paper nếu có thì lại thường và positive nhưng mà nếu mà data làm
[36:12] th nghiệm thì thường là toàn ra negative bạn làm thử nghiệm để thấy để có giáp làm mình không
[36:18] có nhưng là data nó không cân bằng
[36:27] hoặc là giả sử đã predict được buy rồi này thì nó lại là chức năng của nó là như
[36:35] trước mình bảo ví dụ cái bing toket coi như cái ổ khóa chẳng hạn thì cái một cái chất
[36:40] nó gắn vào cái ổ khóa đấy là có thể là nó vặn được cái khóa mà nó chỉ có
[36:46] chui được cái ổ khóa lỗ khóa đây thôi và nó không cho cái chìa khóa chui vào phân loại
[36:51] hai cái đấy cũng là một kiểu n có nhiều bài toán khác nhau nhưng mà cái bài toán cơ
[36:56] bản nhất là một cái target mới và mình cần phải predic một cái chất nào đấy có gắn vào
[37:04] cái target đấy hay không
[37:29] thế ví dụ mà bing thì nó sẽ có những cái chỉ số như nào đính thì mạnh hay yếu
[37:36] nó đo bằng cái ái lực ấ c gì đấy mình không nhớ chv xem mở cái binding cv hoặc
[37:45] là binding dbp ra nó sẽ có nhiều thí nghiệm khác nhau còn nếu mà ức chế về mặt chức
[37:51] năng chẳng hạn thì nó có thể đo bằng cái nồng độ IC50 IC50 là cái nồng độ của cái
[37:56] chất mà nó ức chế một nửa cái chức năng của cái protein đấy
[38:05] ví dụ như rể tu khuẩn chẳng hạn cái đấy thì hơi khác nhưng mà nếu mà muốn chết phuẩn
[38:09] thì nó phải nó phải giết được quá một nửa thì khi đó thì liên tiếp áp dụng như thế
[38:17] thì nó mới chết hết được vi khuẩn
[38:48] thì tuần tiếp thì để làm gì
[39:04] có thể tìm mấy cái bài mới gần đây chỗ tốt đấy xem bài nào mà thấy muốn làm cái
[39:10] gì tương tự thấy thì lại stress vâng nhưng mà những cái bài toán như mình bảo cái bài toán
[39:19] thông dụng nhất đấy là zero short learning
[39:30] tại vì ví dụ một cái target biết rồi nhá có thể tìm được người ta trên Bing DB hoặc
[39:37] là Topcam anh Bình sorry anh Bình em quên mất lịch họp hôm nay em họp ở bên 108 thế
[39:44] là em quên mất lịch họp bây giờ em mới vào sorry anh à không không thì các bạn trình
[39:49] bày xong mấy bài các bạn tìm hiểu rồi đang ngồi nói chuyện vâng ạ vâng thì anh đang bảo
[39:57] mấy bạn có mấy cái bài toán thông dụng thiết kế thuốc xem dạ vâng ạ muốn làm cái gì
[40:04] phải không thì tìm cái bài mới tương tự là mấy cái idea mà cần trước nữa mình bàn thì
[40:11] không biết có bạn nào follow với cả có kết quả gì mới không ạ hôm trước anh có tối
[40:17] trước anh có đưa cái code anh làm cho một cái target cụ thể đấy thì ví dụ như các
[40:23] bạn muốn thì có thể dựa trên cái đấy nhưng mà tu tập các cái target khác thì có thể
[40:31] với các cái target khác nhau thì bao giờ với từng target một mình sẽ có được những cái chất
[40:35] mà họ có kết quả thí nghiệm thí nghiệm thì có sẽ nhiều kiểu khác nhau ví dụ như thúi
[40:40] nghiệm đo nồng độ ức chế thí nghiệm đo bing mạnh hay yếu ờ thì nhưng mà dù gì nữa
[40:47] thì từ các kết quả thí nghiệm đấy cùng một kiểu hạ thì có thể sẽ có vài chất khác
[40:50] nhau thì vài chất đấy có thể mình sẽ phân ra là thành những cái nhóm mà cái chất nó
[40:54] hơi giống nhau thì khi mà trong cùng một nhóm cùng một kiểu thí nghiệm mình sẽ biết là cái
[40:58] nào tôn cái nào thì với từng m cái target đấy mình sẽ thể sinh ra được à một cái
[41:06] ư một cái gì nhỉ một cái tật một cái tật từ chất nọ chất muốn sang chất tốt hơn
[41:15] dạ theo một cái tiêu chí nào đấy dạ thì đấy là một ch đồ thị hả anh không phải
[41:21] dữ liệu thật ấy tức là một target này mà mình đã biết là chất A và chất B nó
[41:29] hơi giống nhau hạ vâng vâng hơi giống nhau ở nhưng mà A và B kết quả thí nghiệm thật
[41:34] nó lại cho một cái chỉ số nào đó ví dụ như Bing chẳng hạn hoặc là ức chế protein
[41:40] vâng dạ đã có số liểu rồi thì mình sẽ biết được là tằng A tốt hơn tằng B thì
[41:46] bây giờ có thể từ đó viết ra cái prom là à vâng từ A thì bây giờ muốn tốt
[41:53] hơn thì như nào thể ra kết quả là B à vâng vâng thế ví dụ từ đấy các bạn
[41:58] ấy có thể build được một cái graph mà có các nốt là các chất còn cạnh là đi từ
[42:04] chất này sang chất kia thì cái gì tốt hơn cái tốt hơn xong rồi dùng cái đấy để cho
[42:08] vào rá cái gì đấy không ạ ờ nếu như coi mỗi một cái protein target là một nút tính
[42:18] chẳng hạn thì mình sẽ có thuận mình có thể thu dượt được rất nhiều cái nút như thế dạ
[42:22] mỗi cái nút đấy thì lại có thể thu thập nhiều cái chắc nhiều compound nhiều molecule vâng có mũi
[42:30] tên đến cái nút đấy và cái cạnh như đấy cạnh đấy nó có thể là cái loại thí nghiệm
[42:35] ví dụ như Vâng đồng độ ức chế và cái property của cái cạnh này nồng độ bằng bao nhiêu
[42:42] thì hoàn toàn có thể xây dựng được những cái Grap như thế vâng vâng nhưng mà giữa các protein
[42:49] với nhau thì giữa các protein với nhau thì có thể dùng một cái beding model gì đấy thì sẽ
[42:56] ra được cái cái độ tương đồng hay là giống nhau vâng thì là như vậy thì có thể là
[43:06] với một protein mới mình chưa có cái chất nào cả nhưng mình có thể xem những cái protein gần
[43:13] giống với nó và xem là cái Vâng nối với cái protein đấy vâng ừ à cái đấy thì cái
[43:23] grap ý tưởng grap của cô thì tuần trước các em cũng về các em tìm hiểu thêm rồi sau
[43:30] đấy thì ở buổi tận trước thì các em bắt đầu thu thập data về và formulate các cái kiểu
[43:37] là bắt đầu xây dựng các cái hướng dataset đấy cô xong rồi Khải tuần trước thì Khải trình bày
[43:43] cho thầy nghe các cái hướng làm set đấy thì thầy có thầy có hướng dẫn Khải lấy co data
[43:50] luôn thì thì Khải ơi em nói qua cho thầy cô nghe cái tiến độ lấy data của bọn c
[43:57] data của bọn em được không về cái cata thì ở hướng thứ nhất ấy là như trước em bảo
[44:04] là cái hướng thứ nhất là cái hướng mà một input ra nhiều output ấy cái dat set đấy thì
[44:10] em có thu được từ ba cái dataset là CBL PEM với lại Papirus thì tổng tất cả dataset của
[44:18] nó là tầm khoảng à khoảng 60.000 sample ạ đấy là hướng thứ nhất ạ còn ở hướng thứ hai
[44:25] là hướng binding dv thì nó sẽ bao gồm một số cái feature khác liên quan đến bing độ bám
[44:32] rồi là độ ức chế thứ như thầy nói vừa nãy thì ở hướng thứ hai thì cái form của
[44:39] nó nó sẽ là một cái input cộng với cái một cái pro target và nó sẽ ra hai cái
[44:45] chất một chất mạnh hơn và một chất yếu hơn thì ở phần này em có được khoảng 280.000 ngì
[44:52] chất để trên ạ và tầm khoảng 30.000 chất để valid cho validation vào test ạ ừ thì đấy là
[45:00] cái thế model của bọn em thì model gì hả Khải ờ hiện tại thì em mới đang đi thu
[45:08] thập data trước ạ tại vì data cái này nó cũng hiếm ạ nó ít ạ nên là em phải
[45:12] em phải lấy từ nhiều nguồn về xong ghép lại ạ ừ ừ ý cô là m à dạ vâng
[45:17] anh Bình nói đi ạ ờ ờ em lấy từ nhiều nguồn lên để xem là chất lượng của nó
[45:21] nhá tại vì thường mình phải lấy những cái mà có thí nghiệm thật đấy vâng em ưu tiên những
[45:26] cái mà nó có thí nghiệm thật ạ ừ cái a cái quên mết cái tên cái gì cái delta
[45:32] sáet thứ ba em bảo ấy nó sẽ có nhiều cái chất lượng khác nhau ấy dạ vâng em cũng
[45:37] có xem thì cái set đấy nó cũng không được như hai cái set đầu tiên là tr ừ ờ
[45:45] model của em input nó sẽ là cái protein cái chất hiện tại và cái prom cần chị sửa sau
[45:52] đấy output của em sẽ là cái chất mà nó thỏa mãn yêu cầu đấy à phải dạ vâng ạ
[45:57] thì input của em sẽ là một cái đầu vào và một cái instruction nó sẽ có ba yêu cầu
[46:02] là tăng cái này giữ nguyên cái này và giảm cái này ví dụ như thế thì đầu ra của
[46:06] em là sẽ là hai cái chất mà tốt hơn ạ thì trước đấy em cũng đã thử một cái
[46:11] set nó đa dạng hơn là hai chất tốt hơn và một chất yếu hơn để model có thể học
[46:15] được tốt hơn nhưng mà nếu mà lấy như thế thì cái số lượng sample nó còn lại rất ít
[46:20] ạ như cái model của em ấy cấu trúc model nó là gì và em tren kiểu gì ờ hiện
[46:26] tại là cái đấy em chưa làm đến ạ thế bây giờ em mới dừng ở cái bức là em
[46:32] craw data thôi còn cái module của em như thế nào và tren như thế nào là chưa có suy
[46:36] nghĩ đúng không dạ vâng em chưa xây dựng cấu trúc model thì các cái model hiện nay người ta
[46:42] thường dùng một có hai hướng một là hai là dùng cái grap diffusion nha cô là kiểu cái grap
[46:50] đầu tiên xong rồi làm cho cái grap đấy nó kiểu nó lung tung đi ấ xong rồi khôi phục
[46:55] lại một cái instruction đấy là nó là condition cho cái diffusion đấy đúng rồi ạ thì nó là hai
[47:02] cái biểu kiểu phổ biến ạ ờ nhưng nó hơi tù mù nhỉ đúng không vâng ờ cái instruction của
[47:09] em là nó có template à tức là chỉ tăng cái này giảm cái kia là nó hoàn toàn có
[47:13] template à vâng nó có một cái temp sẵn ạ nó chỉ thay đổi tên chất và cái số lượng
[47:18] chỉ số thôi ạ tên chất và số lượng chỉ số nếu mà kiểu đấy thì thay vì diffusion mà
[47:24] trên theo kiểu supervive mà xem cái gọi là chia cái instruction đấy thành các cái trường cụ thể và
[47:31] coi như nó là các cái feature input sau đấy thì em chuyên tử s thì nó có chính xác
[47:36] hơn không nhỉ tức một chất coi như là một cái grap thôi xong dùng một cái grap embedding còn
[47:42] các cái instruction kia thì thành tách nó thành các cái feature đầu vào sau đấy thì cho hết nó
[47:47] vào input sau đấy ch theo kiểu ồ encoder decoder mà theo kiểu superv thì nó có chính xác hơn
[47:54] là cái việc là cho vào diffusion không anh thấy mấy cái anh làm anh toàn dùng mấy cái machine
[48:03] classic thôi vâng có khi còn tốt hơn diffusion ấy tại vì là nó ít vâng vâng tại data ít
[48:12] nhưng mà bây giờ thường là các cái bài người ta public thì theo trend đấy cô ta toàn lm
[48:18] với cả diffusion nhưng mà ừ thế thì có sao đâu ví dụ bây giờ em ra một cái em
[48:22] chứng minh là không tốt bằng hay là diffusion tức là generative không tốt bằng những cái phương pháp vai
[48:29] truyền thống thì cũng ok mà miễn là em ra một cái tốt hơn và em chứng minh là nó
[48:34] tốt hơn là được đúng không vâng ạ tại vì LM hay là Zerotip AI cô cảm giác là nó
[48:40] rất là may rủi ấy mình cứ ném vào đấy thôi và nó ra được cái gì mình hoàn toàn
[48:44] không biết ấy bởi vì cái mà nó sinh ra nó chỉ là một cái nó đúng đúng cái phân
[48:48] phối đấy thôi cái phân phối gợi đích đấy thôi còn nó không ra được chính xác cái chất mà
[48:53] mình mong muốn đúng không cái chất mà mình mong muốn nó phải ra một cái là nó nó rất
[48:56] là chính xác ấy còn cái Zen P thì nó chỉ ra một cái mà nó thuộc cái phân phối
[49:01] đấy thôi chứ nó không ra đúng cái instant đấy đúng không nên cô cũng chẳng biết là cái dùng
[49:07] Genip A nó có tốt hơn là cái chuyên cửu vai truyền thống không ấy vâng tại em thầy em
[49:14] thì em nghĩ là cái nếu mà ít data thì cái truyền thống thì nó sẽ có có nó sẽ
[49:21] mạnh hơn đấy hoặc ấy cô nghĩ là nếu mà em dùng generative AI nhiều khi là em dùng thêm
[49:26] cả reinforcement learning ấ maybe maybe là nó work tức là em dùng rất là nhiều cái vòng lặp của
[49:32] generative AI xong rồi cứ chỉnh sửa dần dần làm sao đấy để cho nó đi theo cái hướng ở
[49:38] ở trong cái phân phối đấy đúng đến cái điểm mà mình mong muốn ấy Tiến ạ vâng thì em
[49:46] thấy mấy cái bài mình mấy bài họ làm bây giờ cũng họ cũng đơn giản đấy họ cũng toàn
[49:51] ném vào mấy cái bài benchmark họ toàn ném vào thôi à
[49:58] như bài của năm gần đây thì cái phong cái form chung của nó thường là division hoặc là xong
[50:03] rồi nó sẽ cộng thêm cái RL để tính reward thế à nó cộng thêm RL à ờ thế người
[50:09] ta cũng có cộng thêm RL rồi vâng cũng có một nhiều bài người ta làm thế rồi ạ ừ
[50:13] cô cũng nghĩ là nếu mà dùng generative AI thì phải có cộng thêm RL thì nó mới work còn
[50:18] không thì nó sẽ chả work nhưng mà nếu mà cộng thêm RL ấy thì em sẽ có thêm cái
[50:22] phần mà kiểu là exploration nó như thế nào ấy thì cái chỗ mà exploration thì em có thể add
[50:28] cái huristic của em vào tức là em có thể sinh ra những cái biến thể mới mà nó có
[50:33] khả năng tốt hơn ấy thì cái chỗ đấy em add cái huristic của em vào thì maybe là nó
[50:37] sẽ ra cái tốt hơn đúng không tức là nếu mà em dùng mấy cái generative AI bình thường thì
[50:43] nó cứ đi theo một chiều thôi bây giờ em thêm một cái exploration để nó rẽ ra cái nhanh
[50:47] hướng mới theo những cái tri thức về cái cái cái knowledge về cái cái domain đấy của em nữa
[50:54] ấy để cho nó sinh ra những cái bên cá thể nó khác biệt đi một chút ấy thì cái
[50:58] đấy là cái mình có thể add thêm cái thí thức của mình vào ừ dạ vâng dạ hình như
[51:05] là cái hướng đó là hôm nọ em cũng đã survey có một bài là ở trong người ta đăng
[51:09] trên chước của AI người ta cũng đăng là dùng genetic algorism để có thể đúng rồi ờ cô cũng
[51:16] nghĩ là thế đấy Vâng ạ ừ ừ thì ừ đấy thì cô cũng nghĩ là người ta sẽ làm
[51:21] loanh ngay như thế thì bọn em có thể thử xem những cái method hiện tại xong rồi test với
[51:25] cả các cái bộ raet mà các em collect được ấ xong rồi xem nó bị ở chỗ nào ấy
[51:31] xong các em cải thiện các cái thuật toán RL của họ hoặc là cải thiện các cái process của
[51:37] họ để cho nó tốt hơn đúng không vâng ạ ừ
[51:53] vâng em nghĩ là để cái cũng cứ focus theo cái hướng đấy là ok cho các bạn à ừ
[52:01] thế các bạn là sẽ làm data giống thầy Bình hướng dẫn các thứ sau khi có tập data thì
[52:07] các bạn có thể cài mấy cái baseline đấy xong các bạn test xem là cái nó tốt ở chỗ
[52:11] nào nó không tốt ở chỗ nào để các bạn cải thiện đúng không vâng em nghĩ là cái direction
[52:16] đấy là hợp lý ạ ừ ừ
[52:28] thầy thầy còn câu hỏi gì nữa không thầy
[52:35] ở bên đấy đang nóng hay đang lanh hả anh Bình sắp sờ vào mùa đông đấy ạ bên này
[52:41] nóng lắm anh ạ dạ anh Bình ơi thằng cu con nhà em đợt nó sang chỗ anh ấy ờ
[52:49] vâng nó chuẩn bị đi Mỹ học rồi anh ạ ờ ờ anh cũng nghe nói thế chú mình nhá
[52:55] dạ ôi em cảm ơn mà nó cũng đang thích làm về Grap ấy thế là em bảo thế từ
[53:00] khi nào dọi lên rồi cho giới thiệu xin làm với bác Bình bác Bình cũng làm grap nấu sinh
[53:06] viên làm chứ anh biết vâng bạn ấy cũng làm về heo khe anh ạ thì lúc nào bạn ấy
[53:12] giỏi giỏi thì em xin cho bạn ấy làm cùng với cả các bạn trong nhóm mình giờ thì bạn
[53:17] ấy vẫn đang học thôi bạn ấy giỏi rồi cầu đang còn chập chứng lắm anh ạ à hình như
[53:24] còn được nhận vào cả nữa đúng không ờ vâng đậu cả NUS nữa nhưng mà bạn ấy thích Mỹ
[53:30] hơn là bạn ấy đi Mỹ ạ ừ tị nữa dạ vâng
[53:50] bạn Đăng lần trước sang chỗ Anh ấy thấy Đăng cũng bảo thích sang New Zealand
[53:57] vâng sang năm đăng mới tốt nghiệp à thế à sang năm mới tốt nghiệp vâng sang năm mới tốt
[54:01] nghiệp năm thứ ba thôi sang năm Khải mới Tùng ở đây cũng tốt nghiệp thế thế đấy muốn theo
[54:08] thầy bình thầy cố gắng thầy kiếm tiền
[54:16] bây giờ có rất nhiều học trò cho thầy rồi có học trò áp lực kiếm tiền đồ lên thầy
[54:22] nhiều dạ vâng ạ
[54:30] nay các bạn còn báo cáo gì nữa không Tiến ơi à nãy Tùng có report hai bài cũng làm
[54:38] hai bài benchmark mới ở ICA cho thầy nghe rồi em thấy là format nó cũng tương tự tương tự
[54:46] với bài trước à anh anh anh có anh có gửi vào Zalo mấy cái danh sách tat đang s
[54:56] new là anh phải biết để review thế bảo các bạn ấy xem là có bài nào thấy hay thì
[55:03] bảo anh anh
[55:06] khách về dataset mà nhiều người liên quan vâng
[55:17] máy đ cố gắng build dat nhanh lên xong rồi cài baseline rồi mấy bữa tới report xem là baseline
[55:23] nó có gì hay có gì chưa hay nhá dạ vâng ạ thế cũng cố gắng tìm cái hộ khác
[55:30] nhau ạ đó cũng hơi hơi yếu ạ data hiếm mà mình build được thì nó thành contribution của mình
[55:37] dạ vâng đấy thấy ra data nó cũng khó khá nhiều đấy quan trọng là cái cái lọc nó phải
[55:42] chính xác hoặc là nó phải hợp lý thôi vâng nó cũng nhiều nhưng mà những cái experimental thì nó
[55:47] không nhiều lắm ạ vẫn lại là làm theo yêu cầu đet của mình set phức tạp kiểu như là
[55:51] multi instruction đ thì nó cũng làm về số lượng sample nó cũng giảm đi
[56:11] chỗ đây có ai làm thí nghiệm về hóa sinh nhỉ không anh ạ hôm trước em vừa hỏi mấy
[56:20] bên viện lao phổi trung ương thì cái chị đấy chị bảo chị cũng có làm ở phòng weblast với
[56:26] cả bệnh viện 108 hay là Bạch Mai họ cũng đều có phòng thí nghiệm ấy ạ thì anh cho
[56:32] em biết thông tin cụ thể hơn về cái thí nghiệm anh muốn làm này như thế nào được không
[56:35] ạ của anh thì nó lại ở cái mức thấp hơn mức thấp hơn gồm biên i tức là ví
[56:42] dụ như là có một cái lá nào đấy họ chuyên về một cái protein nào đấy chẳng hạn thì
[56:48] khi mình đưa cho họ một cái chất đấy thì họ họ có thể test là cái chất đấy tác
[56:52] động với cái protein như thế nào ờ đưa cho họ một cái chất thì họ test được là cái
[56:57] chất đấy làm thay đổi protein như thế nào ấy ạ ờ thì tóm lại là mình làm cái gì
[57:01] lại phụ thuộc vào họ cơ chứ không phải là họ phụ thuộc vào mình à tức là họ nghiên
[57:06] c kiếm được một cái gì đấy họ có thể làm được cái thí nghiệm để kiểm chứng cái chất
[57:10] một cái chất bất kỳ đối với cái đấy thì khi đó mình lại sẽ phải collect data thứ mình
[57:17] làm cho họ ờ thế cái chất bất kỳ đấy làm sao để tạo được anh cái chất mới ạ
[57:22] thì sao tạo được ạ à chất nối thì lại cần bên làm hóa họ tạo hoặc là hoặc là
[57:29] mình thuê công ty hóa họ tạo ra hoặc là đơn giản là mình không tạo ra mà mình chỉ
[57:36] à screen những cái thư viện thôi ví dụ như là thư viện có các các cái công ty hóa
[57:41] chất họ sẽ bán sẵn những cái những cái thư viện ví dụ có vài nghìn chất chẳng hạn thì
[57:46] họ sẽ cho mình biết cái cấu trúc mình muốn mua chất nào thì mình sẽ đặt hàng thì họ
[57:50] sẽ gửi cái chất thật đấy về cái để cho test thế cái lap nhận được cái chất đấy thì
[57:56] họ sẽ test các tin chất theo yêu cầu của mình ấy ạ đúng rồi thì khi đó mình sẽ
[58:01] phải dùng mô hình để mình chọn ra cái chất nào để mua để test à em hiểu thế nó
[58:06] ok công ty họ đã có sẵn các cái cấu trúc họ họ bán các chất đấy còn nếu không
[58:12] mình tìm mình tìm ra cấu trúc mới thì mình lại phải thuê người tổng hợp cái chất đấy à
[58:16] thế thì cái phòng kí nghiệp đấy là nó thuộc cái khoa nào của bệnh viện hả anh bệnh viện
[58:22] thì không bệnh viện thì thường là nó ở mức cao hơn rồi họ đã ra được cái chất thành
[58:28] thuốc hoặc là gần thành thuốc để thử nghiệm rồi vâng vâng còn những cái này nó trường đại học
[58:32] thôi trường đại học ở viện nghiên cứu à trường đại học ở viện nghiên cứu biochemistry hoặc là khoa
[58:39] khoản dược ấy em có quen ai ở trường dược hỏi họ làm trường dược à thế em quen nhiều
[58:45] với bệnh viện ấy còn khoa trường dược thì em lại chưa quen thân ai cả thì để em hỏi
[58:49] ạ nhưng mà cái khoa hóa thực phẩm hóa thực phẩm ở Bách Khoa thì anh nghĩ là người ta
[58:56] có phòng láp để làm cái đấy không ạ em cũng thấy người ta có phòng láp thì cũng được
[58:59] thì mình lại phải làm theo cái bài toán của họ vâng thế cái khoa hóa của Bách Khoa thì
[59:05] em khá là quen thân để em hỏi xem họ có hồi trước anh có mấy cái meeting với Trung
[59:14] này với cả Phong này với cả à anh Phong ấy anh Phong ấy ạ vâng thế anh Phong có
[59:21] làm được không ạ thì em ở ngay cạnh phòng của anh Phong nếu mà anh Phong làm được thì
[59:26] là Phong thì làm về kiểu hơi sinh hơi hơi để sinh một chút nhưng là có mấy người làm
[59:32] về
[59:35] sản phẩm chức năng ấy natural product sản phẩm chức năng hồi xưa mấy cái meeting xưa t quên mấy
[59:45] cái tên rồi anh thể lại search lại vâng không thì phải là trường dược đúng không ạ trương vâng
[59:53] anh Trương Quốc Phong thì anh Phong thì em khá là thân thiết thì nếu mà anh Phong làm được
[59:58] thì em nhờ anh Phong easy thực ra thì đáng ra là thỉnh thoảng cũng tinh đấy nhưng mà thấy
[01:00:05] lâu lắm rồi không thấy bên ấy gọi nữa vâng thế là đầu hồi xưa thầy Tháo xuống đây thầy
[01:00:11] về thầy giới thiệu phát vâng một nhóm me tining mỗi lần sau vâng đầu tiên có anh Chu Kỳ
[01:00:18] Sơn này xong rồi anh Chu Kỳ Sơn thì em cũng thân quen ạ anh Chu Kỳ Sơn là hiệu
[01:00:24] trưởng trường hóa được rồi sau đó thì có bạn đấy là Đào Huy Toàn toàn em không biết là
[01:00:35] Biến Thoán vâng em biết anh Sơn với cả anh Phong anh đấy thì em có thể nhờ được thoải
[01:00:42] mái anh Sơn hình như anh cũng làm cái gì đấy đúng không anh Sơn bây giờ anh ấy chỉ
[01:00:49] quản lý thôi chứ anh có làm gì nữa đâu ạ anh ấy là hiệu trưởng anh có làm gì
[01:00:53] nữa đâu còn anh Phong thì vẫn làm ạ anh Phong thì vẫn chết ạ phong Phong Phong là đồng
[01:00:59] đồ nhờ cầu mối mà nên là Vâng nếu mà cần em gặp anh Lênh Phong vâng anh Phong thì
[01:01:06] anh Phong rất là nice anh Phong thì cần support lên hay support lên một thì mình có thể làm
[01:01:11] những cái gì đấy nó có thí nghiệm thật tất nhiên là nó sẽ không những cái jer conference là
[01:01:17] có s nhưng mà nó có thể là những cái gì màạ Vâng dạ nếu cái ngon thì có thể
[01:01:25] đến nature machine intelligence hoặc là Vâng nature dạ vâng hoặc là chắc là phải xem cái bên dược đúng
[01:01:35] không ạ bên trường dược đúng không ạ thì người ta mới làm nhiều đúng không ạ trường dược dược
[01:01:41] hóa cũng được thực ra thì không có biết là công ty thuốc người ta có làm không anh công
[01:01:47] ty thuốc đấy ạ ờ chắc chắn là có nhưng mà vấn đề là làm đấy thì lại phải liên
[01:01:55] quan đến kiểu IP các thứ vâng em đang có cái quan hệ khá là thân thiết với cả bên
[01:02:01] Asraen ạ anh biết công ty không ạ ờ công ty đấy thì nó lập hẳn một cái hub một
[01:02:10] cái eub ở chỗ của em mà em làm với cả asraeneca Việt Nam ấy à thế ở đây là
[01:02:16] cái có lẽ là phù hợp nhất đấy vâng thì với cả họ cũng kết nối em với cả mấy
[01:02:21] cái hub của họ ở trên thế giới ấy tức là bạn Asan đấy nó có một số cái hub
[01:02:25] ở một số nước mà cái hub ở Việt Nam thì nó đặt ở chỗ em nó đặt ở Airfly
[01:02:29] chỗ bọn em thì và em làm với cả bên Asraeneca Việt Nam rất nhiều ấy và họ rất là
[01:02:35] support em thế nhưng mà cái Astraca Việt Nam thì bây giờ họ chỉ bán thuốc là chính thôi chứ
[01:02:40] còn họ chưa có phát triển nghiên cứu về thuốc còn họ có phát triển nghiên cứu về AI thì
[01:02:45] cái phần phát triển nghiên cứu về AI là họ đặt ở chỗ em còn về thuốc là họ bảo
[01:02:49] chưa vì là họ nghĩ là Việt Nam chưa có năng lực để làm cái việc đấy còn họ làm
[01:02:52] cái đấy ở những nước khác ừ vâng phận rồi thì nếu mà cần thì họ dùng luôn bộ phận
[01:02:59] đấy vâng để em hỏi xem là có thể kết nối mình với mấy cái phòng láp về phát triển
[01:03:05] thuốc của họ ở các cái hấp khách trên thế giới được không ấy thì để em hỏi cái đấy
[01:03:09] với mấy bạn ở bên Asra Geneta ok vâng dạ
[01:03:26] vâng ạ dạ vâng chắc mình dừng ở đây ạ vâng em chào thầy cô chào anh Bình chào anh
[01:03:34] Bình chào các em