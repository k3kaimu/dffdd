module dffdd.blockdiagram.model;


enum isModelParameterSet(T) = is(typeof((T t){
    import dffdd.utils.unit;

    // 送信電力が可変であること
    t.txPower = 15.dBm;
    Gain txPower = t.txPower;

    // // 電波暗室内部での動作シミュレーションが可能
    // t.atAnechoicChamber();

    // チャネルとして同軸線路を使用可能
    t.useCoaxialCableAsChannel();

    // INRを設定可能
    t.INR = 10.dB;
}));
